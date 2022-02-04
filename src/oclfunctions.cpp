#include "oclfunctions.h"
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;
static std::string get_error_string(cl_int err)
{
    switch (err) {
        case 0:
            return "CL_SUCCESS";
        case -1:
            return "CL_DEVICE_NOT_FOUND";
        case -2:
            return "CL_DEVICE_NOT_AVAILABLE";
        case -3:
            return "CL_COMPILER_NOT_AVAILABLE";
        case -4:
            return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case -5:
            return "CL_OUT_OF_RESOURCES";
        case -6:
            return "CL_OUT_OF_HOST_MEMORY";
        case -7:
            return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case -8:
            return "CL_MEM_COPY_OVERLAP";
        case -9:
            return "CL_IMAGE_FORMAT_MISMATCH";
        case -10:
            return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case -11:
            return "CL_BUILD_PROGRAM_FAILURE";
        case -12:
            return "CL_MAP_FAILURE";

        case -30:
            return "CL_INVALID_VALUE";
        case -31:
            return "CL_INVALID_DEVICE_TYPE";
        case -32:
            return "CL_INVALID_PLATFORM";
        case -33:
            return "CL_INVALID_DEVICE";
        case -34:
            return "CL_INVALID_CONTEXT";
        case -35:
            return "CL_INVALID_QUEUE_PROPERTIES";
        case -36:
            return "CL_INVALID_COMMAND_QUEUE";
        case -37:
            return "CL_INVALID_HOST_PTR";
        case -38:
            return "CL_INVALID_MEM_OBJECT";
        case -39:
            return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case -40:
            return "CL_INVALID_IMAGE_SIZE";
        case -41:
            return "CL_INVALID_SAMPLER";
        case -42:
            return "CL_INVALID_BINARY";
        case -43:
            return "CL_INVALID_BUILD_OPTIONS";
        case -44:
            return "CL_INVALID_PROGRAM";
        case -45:
            return "CL_INVALID_PROGRAM_EXECUTABLE";
        case -46:
            return "CL_INVALID_KERNEL_NAME";
        case -47:
            return "CL_INVALID_KERNEL_DEFINITION";
        case -48:
            return "CL_INVALID_KERNEL";
        case -49:
            return "CL_INVALID_ARG_INDEX";
        case -50:
            return "CL_INVALID_ARG_VALUE";
        case -51:
            return "CL_INVALID_ARG_SIZE";
        case -52:
            return "CL_INVALID_KERNEL_ARGS";
        case -53:
            return "CL_INVALID_WORK_DIMENSION";
        case -54:
            return "CL_INVALID_WORK_GROUP_SIZE";
        case -55:
            return "CL_INVALID_WORK_ITEM_SIZE";
        case -56:
            return "CL_INVALID_GLOBAL_OFFSET";
        case -57:
            return "CL_INVALID_EVENT_WAIT_LIST";
        case -58:
            return "CL_INVALID_EVENT";
        case -59:
            return "CL_INVALID_OPERATION";
        case -60:
            return "CL_INVALID_GL_OBJECT";
        case -61:
            return "CL_INVALID_BUFFER_SIZE";
        case -62:
            return "CL_INVALID_MIP_LEVEL";
        case -63:
            return "CL_INVALID_GLOBAL_WORK_SIZE";
    }
    return "Unknown OpenCL error";
}

cl_mem oclCreateBuffer(cl_context context, cl_mem_flags flags, size_t size, void* host_ptr)
{
    cl_int err;
    cl_mem mem = clCreateBuffer(context, flags, size, host_ptr, &err);
    if (err != CL_SUCCESS) {
        cerr << "clCreateBuffer Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
    return mem;
}

void oclGetPlatformIDs(cl_uint num_entries, cl_platform_id* platforms, cl_uint* num_platforms)
{
    cl_int err = clGetPlatformIDs(num_entries, platforms, num_platforms);
    if (err != CL_SUCCESS) {
        cerr << "clGetPlatformIDs Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

void oclGetDeviceIDs(cl_platform_id platform,
                     cl_device_type device_type,
                     cl_uint num_entries,
                     cl_device_id* devices,
                     cl_uint* num_devices)
{
    cl_int err = clGetDeviceIDs(platform, device_type, num_entries, devices, num_devices);
    if (err == CL_DEVICE_NOT_FOUND)
        (*num_devices) = 0;
    if (err != CL_SUCCESS && err != CL_DEVICE_NOT_FOUND) {
        cerr << "clGetDeviceIDs Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

cl_context oclCreateContext(const cl_context_properties* properties,
                            cl_uint num_devices,
                            const cl_device_id* devices,
                            void(CL_CALLBACK* pfn_notify)(const char* errinfo,
                                                          const void* private_info,
                                                          size_t cb,
                                                          void* user_data),
                            void* user_data)
{
    cl_int err;
    cl_context context =
      clCreateContext(properties, num_devices, devices, pfn_notify, user_data, &err);
    if (err != CL_SUCCESS) {
        cerr << "clCreateContext Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
    return context;
}

cl_program oclCreateProgramWithSource(cl_context context,
                                      cl_uint count,
                                      const char** strings,
                                      const size_t* lengths)
{
    cl_int err;
    cl_program program = clCreateProgramWithSource(context, count, strings, lengths, &err);
    if (err != CL_SUCCESS) {
        cerr << "clCreateProgramWithSource Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
    return program;
}

void oclBuildProgram(cl_program program,
                     cl_uint num_devices,
                     const cl_device_id* device_list,
                     const char* options,
                     void(CL_CALLBACK* pfn_notify)(cl_program, void* user_data),
                     void* user_data)
{
    cl_int err = clBuildProgram(program, num_devices, device_list, options, pfn_notify, user_data);
    if (err != CL_SUCCESS) {
        cerr << "clBuildProgram Failed: " << get_error_string(err) << endl;
        if (err == CL_BUILD_PROGRAM_FAILURE) {
            for (unsigned int i = 0; i < num_devices; i++) {
                // Determine the size of the log
                size_t log_size;
                clGetProgramBuildInfo(
                  program, device_list[i], CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);

                // Allocate memory for the log
                char* log = (char*)malloc(log_size);

                // Get the log
                clGetProgramBuildInfo(
                  program, device_list[i], CL_PROGRAM_BUILD_LOG, log_size, log, NULL);

                // Print the log
                printf("Log from device %d:\n", i);
                printf("%s\n", log);
                free(log);
            }
        }
        throw runtime_error("");
    }
}

cl_kernel oclCreateKernel(cl_program program, const char* kernel_name)
{
    cl_int err;
    cl_kernel kernel = clCreateKernel(program, kernel_name, &err);
    if (err != CL_SUCCESS) {
        cerr << "clCreateKernel Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
    return kernel;
}

cl_command_queue oclCreateCommandQueue(cl_context context,
                                       cl_device_id device,
                                       cl_command_queue_properties properties)
{
    cl_int err;
    cl_command_queue command_queue = clCreateCommandQueue(context, device, properties, &err);
    if (err != CL_SUCCESS) {
        cerr << "clCreateCommandQueue Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
    return command_queue;
}

void oclGetPlatformInfo(cl_platform_id platform,
                        cl_platform_info param_name,
                        size_t param_value_size,
                        void* param_value,
                        size_t* param_value_size_ret)
{
    cl_int err =
      clGetPlatformInfo(platform, param_name, param_value_size, param_value, param_value_size_ret);
    if (err != CL_SUCCESS) {
        cerr << "clGetDeviceInfo Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

void oclGetDeviceInfo(cl_device_id device,
                      cl_device_info param_name,
                      size_t param_value_size,
                      void* param_value,
                      size_t* param_value_size_ret)
{
    cl_int err =
      clGetDeviceInfo(device, param_name, param_value_size, param_value, param_value_size_ret);
    if (err != CL_SUCCESS) {
        cerr << "clGetDeviceInfo Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

void oclReleaseKernel(cl_kernel kernel)
{
    cl_int err;
    if (kernel != 0) {
        err = clReleaseKernel(kernel);
        if (err != CL_SUCCESS) {
            cerr << "clReleaseKernel Failed: " << get_error_string(err) << endl;
            throw runtime_error("");
        }
    }
}

void oclReleaseCommandQueue(cl_command_queue command_queue)
{
    cl_int err;
    if (command_queue != 0) {
        err = clReleaseCommandQueue(command_queue);
        if (err != CL_SUCCESS) {
            cerr << "clReleaseCommandQueue Failed: " << get_error_string(err) << endl;
            throw runtime_error("");
        }
    }
}

void oclReleaseContext(cl_context context)
{
    cl_int err;
    if (context != 0) {
        err = clReleaseContext(context);
        if (err != CL_SUCCESS) {
            cerr << "clReleaseContext Failed: " << get_error_string(err) << endl;
            throw runtime_error("");
        }
    }
}

void oclEnqueueWriteBuffer(cl_command_queue command_queue,
                           cl_mem buffer,
                           cl_bool blocking_write,
                           size_t offset,
                           size_t cb,
                           const void* ptr,
                           cl_uint num_events_in_wait_list,
                           const cl_event* event_wait_list,
                           cl_event* event)
{
    cl_int err = clEnqueueWriteBuffer(command_queue,
                                      buffer,
                                      blocking_write,
                                      offset,
                                      cb,
                                      ptr,
                                      num_events_in_wait_list,
                                      event_wait_list,
                                      event);
    if (err != CL_SUCCESS) {
        cerr << "clEnqueueWriteBuffer Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

void oclFinish(cl_command_queue command_queue)
{
    cl_int err = clFinish(command_queue);
    if (err != CL_SUCCESS) {
        cerr << "clFinish Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

void oclSetKernelArg(cl_kernel kernel, cl_uint arg_index, size_t arg_size, const void* arg_value)
{
    cl_int err = clSetKernelArg(kernel, arg_index, arg_size, arg_value);
    if (err != CL_SUCCESS) {
        cerr << "clSetKernelArg Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

void oclEnqueueNDRangeKernel(cl_command_queue command_queue,
                             cl_kernel kernel,
                             cl_uint work_dim,
                             const size_t* global_work_offset,
                             const size_t* global_work_size,
                             const size_t* local_work_size,
                             cl_uint num_events_in_wait_list,
                             const cl_event* event_wait_list,
                             cl_event* event)
{
    cl_int err = clEnqueueNDRangeKernel(command_queue,
                                        kernel,
                                        work_dim,
                                        global_work_offset,
                                        global_work_size,
                                        local_work_size,
                                        num_events_in_wait_list,
                                        event_wait_list,
                                        event);
    if (err != CL_SUCCESS) {
        cerr << "clEnqueueNDRangeKernel Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

void oclEnqueueReadBuffer(cl_command_queue command_queue,
                          cl_mem buffer,
                          cl_bool blocking_read,
                          size_t offset,
                          size_t cb,
                          void* ptr,
                          cl_uint num_events_in_wait_list,
                          const cl_event* event_wait_list,
                          cl_event* event)
{
    cl_int err = clEnqueueReadBuffer(command_queue,
                                     buffer,
                                     blocking_read,
                                     offset,
                                     cb,
                                     ptr,
                                     num_events_in_wait_list,
                                     event_wait_list,
                                     event);
    if (err != CL_SUCCESS) {
        cerr << "clEnqueueReadBuffer Failed: " << get_error_string(err) << endl;
        throw runtime_error("");
    }
}

void clearbuf(cl_mem buf)
{
    cl_int err;
    if (buf != 0) {
        err = clReleaseMemObject(buf);
        if (err != CL_SUCCESS) {
            cerr << "clReleaseMemObject Failed: " << get_error_string(err) << endl;
            throw runtime_error("");
        }
    }
}
