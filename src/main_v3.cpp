#include "RangeIterator.h"
#include "async_search.h"
#include "bidirect_search.h"
#include "bit4ops.h"
#include "blockify.h"
#include "bulge_logic.h"
#include "ceildiv.h"
#include "parse_input.h"
#include "search.h"
#include "usage.h"
#include <cassert>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

#define CAS_OFFINDER_VERSION "3.1.0"

static DeviceType get_dev_ty(char c)
{
    switch (c) {
        case 'C': return CPU;
        case 'G': return GPU;
        case 'A': return ACCEL;
    }
    throw runtime_error(
      "not a valid device name, must be one of 'G', 'C', 'A'");
}
struct OutDatav3
{
    ostream* out_str_ptr;
    InFileInfo input;
    size_t aug_size;
    std::vector<bulge_pair> augmented_patterns;
    string mirrored_all_cmps;
    string clipped_pattern;
    string rev_clipped_pattern;
    bool reversed_pam;
    std::mutex lock;
};

static void async_callback(const GenomeMatch* gm, void* user_data)
{
    OutDatav3* data = (OutDatav3*)(user_data);
    ostream& outs = *data->out_str_ptr;
    size_t n_augs =
      data->augmented_patterns.size() / (data->input.num_patterns * 2);
    bool is_mirrored = (gm->pattern_idx / n_augs) % 2;
    size_t orig_pidx = (gm->pattern_idx / n_augs) / 2;
    std::string dna(gm->dna_match);
    std::string rna(data->augmented_patterns.at(gm->pattern_idx).first);
    bulge_augment augment = data->augmented_patterns.at(gm->pattern_idx).second;
    bulge_info bi = get_bulge_info(dna,
                                   rna,
                                   augment,
                                   gm->chrom_loc,
                                   data->input.dna_bulges);

    size_t psize = data->clipped_pattern.size();
    if ((!is_mirrored &&
         ((data->reversed_pam &&
           !matches(bi.dna.data(), data->clipped_pattern.data(), psize)) ||
          (!data->reversed_pam &&
           !matches(bi.dna.data() + bi.dna.size() - psize,
                    data->clipped_pattern.data(),
                    psize)))) ||
        (is_mirrored &&
         ((data->reversed_pam && !matches(bi.dna.data() + bi.dna.size() - psize,
                                          data->rev_clipped_pattern.data(),
                                          psize)) ||
          (!data->reversed_pam && !matches(bi.dna.data(),
                                           data->rev_clipped_pattern.data(),
                                           psize))))) {
        return;
    }
    string id = to_string(orig_pidx);
    if (data->input.ids) {
        id = string(data->input.ids[orig_pidx]);
    }
    if (is_mirrored) {
        bi.dna = reverse_compliment(bi.dna);
        bi.rna = reverse_compliment(bi.rna);
    }
    char dir = is_mirrored ? '-' : '+';
    indicate_mismatches_dna(bi.dna, bi.rna);
    stringstream ss;
    ss << id << '\t' << get_bulge_type_name(augment.bulge_type) << '\t'
       << bi.rna << '\t' << bi.dna << '\t' << gm->chrom_name << '\t' << bi.loc
       << '\t' << dir << '\t' << gm->mismatches << '\t' << augment.bulge_size
       << '\n';
    data->lock.lock();
    outs << ss.str();
    data->lock.unlock();
}
static const char usage_str[] = \
        "Cas-OFFinder " CAS_OFFINDER_VERSION "\n"\
        "\n"\
        "Copyright (c) 2021 Jeongbin Park and Sangsu Bae"  "\n"\
        "Website: "  CAS_OFFINDER_HOMEPAGE_URL  "\n"\
        "\n"\
        "Usage: cas-offinder [options] {input_filename|-} "\
          "{C|G|A}[device_id(s)] {output_filename|-}"\
        "\n"\
        "(C: using CPUs, G: using GPUs, A: using accelerators)"  "\n"\
        "\n"\
        "Options"  "\n"\
        "  --summary <file>        Print summary table to the specified file."\
        "\n"\
        "\n"\
        "Example input file (DNA bulge 2, RNA bulge 1):"  "\n"\
        "/var/chromosomes/human_grch38"  "\n"\
        "NNNNNNNNNNNNNNNNNNNNNRG 2 1"  "\n"\
        "GGCCGACCTGTCGCTGACGCNNN 5"  "\n"\
        "CGCCAGCGTCAGCGACAGGTNNN 5"  "\n"\
        "ACGGCGCCAGCGTCAGCGACNNN 5"  "\n"\
        "GTCGCTGACGCTGGCGCCGTNNN 5"  "\n"\
        "\n"\
        "Available device list:"  "\n";

int main(int argc, char** argv)
{
    if (argc != 4) {
        cerr << usage_str;
        SearchFactory* fact = create_search_factory(ANY_DEV);
        print_device_information(fact);
        return 1;
    }
    if (argc != 4) {
        cerr << "requires 3 CLI arguments, in_file, device, out_file\n";
        return 1;
    }
    const char* in_fname = argv[1];
    if (strlen(argv[2]) != 1) {
        cerr << "device should be 'C','G', or 'A'\n";
        return 1;
    }
    char device_chr = argv[2][0];
    const char* out_fname = argv[3];
    bool is_reversed_pam = false;

    auto start = std::chrono::system_clock::now();

    InFileInfo input = read_file(in_fname);
    DeviceType device_ty = get_dev_ty(device_chr);
    ofstream out_file;
    vector<string> orig_cmps(input.num_patterns);
    for (size_t i : range(input.num_patterns)) {
        orig_cmps[i] = string(input.compares + i * input.pattern_size,
                              input.compares + (i + 1) * input.pattern_size);
    }
    vector<string> mir_orig_cmps(input.num_patterns * 2);
    for (size_t i : range(input.num_patterns)) {
        mir_orig_cmps[i * 2] = orig_cmps[i];
        mir_orig_cmps[i * 2 + 1] = reverse_compliment(orig_cmps[i]);
    }
    string orig_pattern(input.pattern, input.pattern + input.pattern_size);
    string clipped_pattern = orig_pattern;
    if (input.dna_bulges != 0 || input.rna_bulges != 0) {
        if (orig_pattern.front() != 'N' && orig_pattern.back() != 'N') {
            cerr << "Pattern must begin or end with 'N', will affect pattern "
                    "application. \n";
            return 1;
        }
        if (orig_pattern.find_first_not_of('N') != string::npos) {
            if (orig_pattern.front() != 'N') {
                clipped_pattern =
                  orig_pattern.substr(0, orig_pattern.find_first_not_of('N'));
                is_reversed_pam = true;
            } else {
                clipped_pattern =
                  orig_pattern.substr(orig_pattern.find_last_not_of('N'));
                is_reversed_pam = false;
            }
        } else {
            clipped_pattern = "";
        }
    }
    std::vector<bulge_pair> augmented_patterns = augment_patterns_with_bulges(
      mir_orig_cmps, input.dna_bulges, input.rna_bulges);
    string augmented_cmp;
    for (bulge_pair bp : augmented_patterns) {
        augmented_cmp += bp.first;
    }
    size_t aug_size = augmented_patterns.at(0).first.size();

    ostream* out_str_ptr;
    if (strlen(out_fname) == 1 && out_fname[0] == '-') {
        out_str_ptr = &cout;
    } else {
        out_file.open(out_fname);
        out_str_ptr = &out_file;
    }
    (*out_str_ptr) << "##Generated by Cas-OFFinder 3.1.0\n";
    (*out_str_ptr) << "#Id\tBulge "
                      "Type\tcrRNA\tDNA\tChromosome\tLocation\tDirection\tMisma"
                      "tches\tBulge Size"
                   << endl;
    OutDatav3 out_data{
        .out_str_ptr = out_str_ptr,
        .input = input,
        .aug_size = aug_size,
        .augmented_patterns = augmented_patterns,
        .mirrored_all_cmps = augmented_cmp,
        .clipped_pattern = clipped_pattern,
        .rev_clipped_pattern = reverse_compliment(clipped_pattern),
        .reversed_pam = is_reversed_pam,
                .lock=std::mutex()
    };
    async_search(input.genome_path,
                 device_ty,
                 augmented_cmp.data(),
                 aug_size,
                 augmented_patterns.size(),
                 input.mismatches,
                 nullptr,
                 &out_data,
                 async_callback);

    auto end = std::chrono::system_clock::now();
    auto duration = (end - start);
    auto micros =
      std::chrono::duration_cast<std::chrono::microseconds>(duration);
    cerr << "successfully finished in: " << micros.count() / 1000000.
         << " seconds\n";
    return 0;
}
