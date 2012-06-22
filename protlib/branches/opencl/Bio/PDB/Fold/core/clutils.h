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
 * This file contains the definition of the functions that control OpenCL
 * computations.
 *
 */

#ifndef __CLUTILS_H__
#define __CLUTILS_H__

#include "datatypes.h"
#include "errorcodes.h"
#include <CL/cl.h>

/**
 * Intializes OpenCL context, devices and command queue.
 * @param clcontext cl_context*
 * @param cldevices cl_device_id**
 * @param clqueue cl_command_queue*
 * @return ErrorCode
 */
ErrorCode init_opencl(cl_context *clcontext, cl_device_id **cldevices, cl_command_queue *clqueue);

/**
 * Frees coords, distance and distance mask memory objects and devices list.
 * @param cldevices cl_device_id**
 * @param clcoords cl_mem*
 * @param cldistmask cl_mem*
 * @param cldistances cl_mem*
 * @return ErrorCode
 */
ErrorCode finish_opencl(cl_device_id **cldevices, cl_mem *clcoords, cl_mem *cldistmask, cl_mem *cldistances);

/**
 * Frees coords, distance and distance mask memory objects and devices list.
 * @param clcoords cl_mem*
 * @param cldistmask cl_mem*
 * @param cldistances cl_mem*
 * @param clcontext cl_context
 * @param coords FloatValue*
 * @param distmask BoolValue*
 * @param distances FloatValue*
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode set_opencl_memory(cl_mem *clcoords, cl_mem *cldistmask, cl_mem *cldistances, cl_context clcontext, FloatValue *coords, BoolValue *distmask, FloatValue *distances, IndexValue natoms);


/**
 * Frees coords, distance and distance mask memory objects and devices list.
 * @param clprogram cl_program*
 * @param cldistkernel cl_kernel*
 * @param cldistsqkernel cl_kernel*
 * @param clinvdistkernel cl_kernel*
 * @param clcontext cl_context
 * @return ErrorCode
 */
ErrorCode compile_opencl_kernels(cl_program *clprogram, cl_kernel *cldistkernel, cl_kernel *cldistsqkernel, cl_kernel *clinvdistkernel, cl_context clcontext);

// Runs the distkernel onthe provided data.
ErrorCode run_opencl_distkernel(cl_command_queue clqueue, cl_kernel cldistkernel, cl_mem clcoords, cl_mem cldistmask, cl_mem cldistances, IndexValue natoms, FloatValue *distances);

#endif

