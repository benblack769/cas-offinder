#pragma once
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <iostream>

using namespace std;

int read_fasta(string &filepath, vector<string> &chrnames, string &content, vector<uint64_t> &chrpos);