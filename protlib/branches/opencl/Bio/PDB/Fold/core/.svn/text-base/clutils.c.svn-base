/* *
 * The C Protein Folding Library.
 * Copyright (C) 2009 Andres Colubri.
 * Contact: andres.colubri 'AT' gmail.com
 *
 * This library was written at the Institute of Biophysical Dynamics at the University of Chicago.
 * Gordon Center for Integrated Science, W101. 929 East 57th Street, Chicago, IL 60637, USA.
 * Homepage: http://ibd.uchicago.edu/
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

/**
 * This file contains the implementation of the functions that control OpenCL
 * computations.
 *
 */

#include <stdio.h>
#include "clutils.h"

/*
static const int src_length = 5;
static const char* dist_calc_src[] = {
  "__kernel void distsq_calc(__global float4 *coords, __global float *distances)\n",
  "{\n",
  "unsigned int i = get_global_id(0);\n",
  "unsigned int j = get_global_id(1);\n",
  "}\n"
};
*/

#define SIZE 2048 // number of elements in the vectors to be added

cl_mem GPUOutputVector;

// OpenCL source code
const char* OpenCLSource[] = {
  "__kernel void VectorAdd(__global int* c, __global int* a,\n",
  " __global int* b, __global float4 *coords, __global char *distmask, __global float *distances, uint natoms)\n",
  "{\n",
  " // Index of the elements to add\n",
  " unsigned int n = get_global_id(0);\n",
  " unsigned int m = get_global_id(1);\n"
  " // Sum the nth element of vectors a and b and store in c\n",
  " c[n] = natoms;\n",
  "}\n"
};

// Some interesting data for the vectors
int InitialData1[14] = {37,50,54,50,56,12,37,45,77,81,92,56,-22,-4};
int InitialData2[14] = {35,51,54,58,55,32,-5,42,34,33,16,44,55,14};

static const int src_length = 44;
static const char* dist_calc_src[] = {
"__kernel void distsq_calc(__global float4 *coords, __global char *distmask, __global float *distances, uint natoms)\n",
"{\n",
"    unsigned int i = get_global_id(0);\n",
"    unsigned int j = get_global_id(1);\n",
"\n",
"    float4 vi, vj;\n",
"\n",
"    vi = coords[i];\n",
"    vj = coords[j];\n",
"\n",
"    int idx = i * natoms + j;\n",
"\n",
"    distances[idx] = distmask[idx] * dot(vi, vj);\n",
"}\n",
"\n",
"__kernel void invdist_calc(__global float4 *coords, __global char *distmask, __global float *distances, uint natoms)\n",
"{\n",
"    unsigned int i = get_global_id(0);\n",
"    unsigned int j = get_global_id(1);\n",
"\n",
"    float4 vi, vj;\n",
"\n",
"    vi = coords[i];\n",
"    vj = coords[j];\n",
"\n",
"    int idx = i * natoms + j;\n",
"\n",
"    distances[idx] = distmask[idx] * 1.0 / distance(vi, vj);\n",
"}\n",
"\n",
"__kernel void dist_calc(__global float *outdist, __global float4 *coords, __global char *distmask, __global float *distances, uint natoms)\n",
"{\n",
"    unsigned int i = get_global_id(0);\n",
"    unsigned int j = get_global_id(1);\n",
"\n",
"    float4 vi, vj;\n",
"\n",
"    vi = coords[i];\n",
"    vj = coords[j];\n",
"\n",
"    int idx = i * natoms + j;\n",
"\n",
"    distances[idx] = distance(vi, vj);\n",
"}\n"
};

ErrorCode init_opencl(cl_context *clcontext, cl_device_id **cldevices, cl_command_queue *clqueue)
{
    // Create a compute context with GPU device:
    *clcontext = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);

    // Then we can get the list of GPU devices associated with this context:
    size_t databytes;
    clGetContextInfo(*clcontext, CL_CONTEXT_DEVICES, 0, NULL, &databytes);
    *cldevices = (cl_device_id*)malloc(databytes);
    clGetContextInfo(*clcontext, CL_CONTEXT_DEVICES, databytes, *cldevices, NULL);

    // Create a work-queue on the first OpenCL device:
    *clqueue = clCreateCommandQueue(*clcontext, (*cldevices)[0], 0, NULL);

    printf("Initializing OPENCL...\n");

    return NO_ERROR;
}


ErrorCode finish_opencl(cl_device_id **cldevices, cl_mem *clcoords, cl_mem *cldistmask, cl_mem *cldistances)
{
    printf("Finalizing OPENCL...\n");

    if (*clcoords != 0) clReleaseMemObject(*clcoords);
    if (*cldistmask != 0) clReleaseMemObject(*cldistmask);
    if (*cldistances != 0) clReleaseMemObject(*cldistances);

    clReleaseMemObject(GPUOutputVector);

    if (*cldevices != NULL) free(*cldevices);

    return NO_ERROR;
}

ErrorCode set_opencl_memory(cl_mem *clcoords, cl_mem *cldistmask, cl_mem *cldistances, cl_context clcontext, FloatValue *coords, BoolValue *distmask, FloatValue *distances, IndexValue natoms)
{
    // Create the IO variables on the compute device, mapping them to the coords,
    // distmask and distances arrays on CPU memory (because of CL_MEM_USE_HOST_PTR).
    *clcoords = clCreateBuffer(clcontext, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(FloatValue) * 4 * natoms, coords, NULL);
    *cldistmask = clCreateBuffer(clcontext, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, sizeof(BoolValue)* natoms * natoms, distmask, NULL);
    *cldistances = clCreateBuffer(clcontext, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, sizeof(FloatValue) * natoms * natoms, distances, NULL);

    // Allocate output memory on GPU
    GPUOutputVector = clCreateBuffer(clcontext, CL_MEM_WRITE_ONLY, sizeof(FloatValue) * natoms * natoms, NULL, NULL);


    printf("Creating OpenCL device variables: %d %d %d\n", *clcoords, *cldistmask, *cldistances);

    return NO_ERROR;
}

ErrorCode compile_opencl_kernels(cl_program *clprogram, cl_kernel *cldistkernel, cl_kernel *cldistsqkernel, cl_kernel *clinvdistkernel, cl_context clcontext)
{
    cl_int clerror;

    // Create the program object.
    *clprogram = clCreateProgramWithSource(clcontext, src_length, dist_calc_src, NULL, &clerror);
    if (clerror != CL_SUCCESS)
    {
        printf("Error during creation of OpenCL program.\n");
    }
    else printf("OpenCL program created succesfully.\n");

    // Build the compute program executable.
    clerror = clBuildProgram(*clprogram, 0, NULL, NULL, NULL, NULL);
    if (clerror != CL_SUCCESS)
    {
        printf("Error during building of OpenCL program.\n");
    }
    else printf("OpenCL program built succesfully.\n");


    // Create the compute kernels:
    *cldistkernel = clCreateKernel(*clprogram, "dist_calc", &clerror);
    if (clerror != CL_SUCCESS)
    {
        printf("Error during creation of dist_calc OpenCL kernel.\n");
    }
    else printf("dist_calc OpenCL kernel created succesfully: %d\n", *cldistkernel);

    *cldistsqkernel = clCreateKernel(*clprogram, "distsq_calc", &clerror);
    if (clerror != CL_SUCCESS)
    {
        printf("Error during creation of distsq_calc OpenCL kernel.\n");
    }
    else printf("distsq_calc OpenCL kernel created succesfully: %d\n", *cldistsqkernel);

    *clinvdistkernel = clCreateKernel(*clprogram, "invdist_calc", &clerror);
    if (clerror != CL_SUCCESS)
    {
        printf("Error during creation of invdist_calc OpenCL kernel.\n");
    }
    else printf("invdist_calc OpenCL kernel created succesfully: %d\n", *clinvdistkernel);

    return NO_ERROR;
}


ErrorCode run_opencl_distkernel(cl_command_queue clqueue, cl_kernel cldistkernel, cl_mem clcoords, cl_mem cldistmask, cl_mem cldistances, IndexValue natoms, FloatValue *distances)
{
    cl_int clerror;
    size_t clrange[2] = {natoms, natoms};

    cl_uint n = natoms;
    clerror = clSetKernelArg(cldistkernel, 0, sizeof(cl_mem), (void*)&GPUOutputVector);
    clerror |= clSetKernelArg(cldistkernel, 1, sizeof(cl_mem), (void *)&clcoords);
    clerror |= clSetKernelArg(cldistkernel, 2, sizeof(cl_mem), (void *)&cldistmask);
    clerror |= clSetKernelArg(cldistkernel, 3, sizeof(cl_mem), (void *)&cldistances);
    clerror |= clSetKernelArg(cldistkernel, 4, sizeof(cl_uint), (void*)&n);
    if (clerror != CL_SUCCESS)
    {
        printf("Error while setting arguments for distkernel.\n");
    }
    else printf("OpenCL arguments OK!\n");


    // When the numBodies / thread block size is < # multiprocessors 
    // (16 on G80), the GPU is underutilized. For example, with 256 threads per
    // block and 1024 bodies, there will only be 4 thread blocks, so the 
    // GPU will only be 25% utilized.  To improve this, we use multiple threads
    // per body.  We still can use blocks of 256 threads, but they are arranged
    // in q rows of p threads each.  Each thread processes 1/q of the forces 
    // that affect each body, and then 1/q of the threads (those with 
    // threadIdx.y==0) add up the partial sums from the other threads for that 
    // body.  To enable this, use the "--p=" and "--q=" command line options to
    // this example.  e.g.: "nbody.exe --n=1024 --p=64 --q=4" will use 4 
    // threads per body and 256 threads per block. There will be n/p = 16 
    // blocks, so a G80 GPU will be 100% utilized.
    size_t global_work_size[2];
    size_t local_work_size[2];
    // set work-item dimensions
    local_work_size[0] = 4;
    local_work_size[1] = 4;
    global_work_size[0]= natoms;
    global_work_size[1]= natoms;

    printf("Running distance on the GPU...\n");
    clerror = clEnqueueNDRangeKernel(clqueue, cldistkernel, 2, NULL, clrange, NULL, 0, NULL, NULL);
    if (clerror != CL_SUCCESS)
    {
        if (clerror == CL_INVALID_PROGRAM_EXECUTABLE) printf("CL_INVALID_PROGRAM_EXECUTABLE\n");
        else if (clerror == CL_INVALID_COMMAND_QUEUE) printf("CL_INVALID_COMMAND_QUEUE\n");
        else if (clerror == CL_INVALID_KERNEL) printf("CL_INVALID_KERNEL\n");
        else if (clerror == CL_INVALID_CONTEXT) printf("CL_INVALID_CONTEXT\n");
        else if (clerror == CL_INVALID_KERNEL_ARGS) printf("CL_INVALID_KERNEL_ARGS\n");
        else if (clerror == CL_INVALID_WORK_DIMENSION) printf("CL_INVALID_WORK_DIMENSION\n");
        else if (clerror == CL_INVALID_WORK_GROUP_SIZE) printf("CL_INVALID_WORK_GROUP_SIZE\n");
        else if (clerror == CL_INVALID_WORK_ITEM_SIZE) printf("CL_INVALID_WORK_ITEM_SIZE\n");
        else if (clerror == CL_INVALID_GLOBAL_OFFSET) printf("CL_INVALID_GLOBAL_OFFSET\n");
        else if (clerror == CL_OUT_OF_RESOURCES) printf("CL_OUT_OF_RESOURCES\n");
        else if (clerror == CL_MEM_OBJECT_ALLOCATION_FAILURE) printf("CL_MEM_OBJECT_ALLOCATION_FAILURE\n");
        else if (clerror == CL_INVALID_EVENT_WAIT_LIST) printf("CL_INVALID_EVENT_WAIT_LIST\n");
        else if (clerror == CL_OUT_OF_HOST_MEMORY) printf("CL_OUT_OF_HOST_MEMORY\n");
        else printf("UNKNOWN CL ERROR: %d\n", clerror);


    }
    else 
    {
       printf("Excecution succesful!\n");       

       FloatValue *HostOutputVector = NULL;
       create_distances_array(&HostOutputVector, natoms);
  
       // Copy the output in GPU memory back to CPU memory
       clEnqueueReadBuffer(clqueue, cldistances, CL_TRUE, 0, sizeof(FloatValue) * natoms * natoms, distances, 0, NULL, NULL);

       int i;
       for (i = 0; i < 10; i++)
           printf("%i %f\n", i, HostOutputVector[i]);
       printf("\n");
  
       delete_distances_array(&HostOutputVector, natoms);
    }




/*
  // Here's two source vectors in CPU (Host) memory
  int HostVector1[SIZE];
  int HostVector2[SIZE];

  // Lets initialize them with some interesting repeating data
  int c;
  for(c = 0;c<SIZE;c++) {
    HostVector1[c] = InitialData1[c%14];
    HostVector2[c] = InitialData2[c%14];
  }

  // We want to run our OpenCL on our CUDA-enable NVIDIA hardware,
  // so we create a GPU type context
  cl_context GPUContext = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, NULL);

  // Then we can get the list of GPU devices associated with this context
  size_t ParmDataBytes;
  clGetContextInfo(GPUContext, CL_CONTEXT_DEVICES, 0, NULL, &ParmDataBytes);
  cl_device_id* GPUDevices = (cl_device_id*)malloc(ParmDataBytes);
  clGetContextInfo(GPUContext, CL_CONTEXT_DEVICES, ParmDataBytes, GPUDevices, NULL);

  // And create a command-queue on the first GPU device
  cl_command_queue GPUCommandQueue = clCreateCommandQueue(GPUContext, GPUDevices[0], 0, NULL);

  // Allocate GPU memory for the source vectors and initialize with the CPU
  // memory
  cl_mem GPUVector1 = clCreateBuffer(GPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * SIZE, HostVector1, NULL);
  cl_mem GPUVector2 = clCreateBuffer(GPUContext, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(int) * SIZE, HostVector2, NULL);

  // Allocate output memory on GPU
  cl_mem GPUOutputVector;
  GPUOutputVector = clCreateBuffer(GPUContext, CL_MEM_WRITE_ONLY, sizeof(int) * SIZE, NULL, NULL);

  // We are now ready to run the kernel code on the GPU
  // OpenCL supports runtime code compilation. First we build the OpenCL
  // program from source code
  cl_program OpenCLProgram = clCreateProgramWithSource(GPUContext, 8, OpenCLSource, NULL, NULL);
  clBuildProgram(OpenCLProgram,0,NULL,NULL,NULL,NULL);

  // Then we can create a handle to the compiled OpenCL function (Kernel)
  cl_kernel OpenCLVectorAdd = clCreateKernel(OpenCLProgram, "VectorAdd", NULL);
  
  // In the next step we associate the GPU memory with the Kernel arguments
  cl_uint n = natoms;
  clSetKernelArg(OpenCLVectorAdd, 0, sizeof(cl_mem), (void*)&GPUOutputVector);
  clSetKernelArg(OpenCLVectorAdd, 1, sizeof(cl_mem), (void*)&GPUVector1);
  clSetKernelArg(OpenCLVectorAdd, 2, sizeof(cl_mem), (void*)&GPUVector2);
  clSetKernelArg(OpenCLVectorAdd, 3, sizeof(cl_mem), (void*)&clcoords);
  clSetKernelArg(OpenCLVectorAdd, 4, sizeof(cl_mem), (void*)&cldistmask);
  clSetKernelArg(OpenCLVectorAdd, 5, sizeof(cl_mem), (void*)&cldistances);
  clSetKernelArg(OpenCLVectorAdd, 6, sizeof(cl_uint), (void*)&n);

  // Then we launch the Kernel on the GPU
  size_t WorkSize[1] = {SIZE, SIZE};
  cl_int clerror = clEnqueueNDRangeKernel(GPUCommandQueue, OpenCLVectorAdd, 2, NULL, WorkSize, NULL, 0, NULL, NULL);
    if (clerror != CL_SUCCESS)
    {
        if (clerror == CL_INVALID_PROGRAM_EXECUTABLE) printf("CL_INVALID_PROGRAM_EXECUTABLE\n");
        else if (clerror == CL_INVALID_COMMAND_QUEUE) printf("CL_INVALID_COMMAND_QUEUE\n");
        else if (clerror == CL_INVALID_KERNEL) printf("CL_INVALID_KERNEL\n");
        else if (clerror == CL_INVALID_CONTEXT) printf("CL_INVALID_CONTEXT\n");
        else if (clerror == CL_INVALID_KERNEL_ARGS) printf("CL_INVALID_KERNEL_ARGS\n");
        else if (clerror == CL_INVALID_WORK_DIMENSION) printf("CL_INVALID_WORK_DIMENSION\n");
        else if (clerror == CL_INVALID_WORK_GROUP_SIZE) printf("CL_INVALID_WORK_GROUP_SIZE\n");
        else if (clerror == CL_INVALID_WORK_ITEM_SIZE) printf("CL_INVALID_WORK_ITEM_SIZE\n");
        else if (clerror == CL_INVALID_GLOBAL_OFFSET) printf("CL_INVALID_GLOBAL_OFFSET\n");
        else if (clerror == CL_OUT_OF_RESOURCES) printf("CL_OUT_OF_RESOURCES\n");
        else if (clerror == CL_MEM_OBJECT_ALLOCATION_FAILURE) printf("CL_MEM_OBJECT_ALLOCATION_FAILURE\n");
        else if (clerror == CL_INVALID_EVENT_WAIT_LIST) printf("CL_INVALID_EVENT_WAIT_LIST\n");
        else if (clerror == CL_OUT_OF_HOST_MEMORY) printf("CL_OUT_OF_HOST_MEMORY\n");
        else printf("UNKNOWN CL ERROR: %d\n", clerror);


    }
    else printf("Excecution succesful!\n");

  // Here’s a vector in CPU memory the output
  int HostOutputVector[SIZE];
  
  // Copy the output in GPU memory back to CPU memory
  clEnqueueReadBuffer(GPUCommandQueue, GPUOutputVector, CL_TRUE, 0, SIZE*sizeof(int), HostOutputVector, 0, NULL, NULL);

  // We are finished with GPU memory so we can free it
  clReleaseMemObject(GPUVector1);
  clReleaseMemObject(GPUVector2);
  clReleaseMemObject(GPUOutputVector);
  free(GPUDevices);

  // Print out the results for fun.
  // We are simply casting the numeric result to a char
  // and printing it to the console
  for(c = 0; c < 305;c++)
    printf("%d|",HostOutputVector[c]);
*/

    return NO_ERROR;
}




/*
// run
// Set the arguments for the kernel

    // Initialize NDRange
    clrange[0] = natoms;
    clrange[1] = natoms;

//size_t WorkSize[1] = {SIZE};
//clEnqueueNDRangeKernel(GPUCommandQueue, OpenCLVectorAdd,

// Then we launch the Kernel on the CL device.

// execute kernel
//clExecuteKernel(queue, kernel, NULL, range, NULL, 0, NULL);
*/


