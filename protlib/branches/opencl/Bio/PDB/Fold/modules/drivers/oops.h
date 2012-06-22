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
 * Default simulation driver.
 *
 */

#ifndef __OOPS_H__
#define __OOPS_H__

#include "../../core/datatypes.h"

#ifdef __USE_OPENCL_CODEPATH__
#include "../../core/clutils.h"
#endif

ErrorCode set_system_for_oops(struct AtomData *_atoms, FloatValue *_coords, 
                                                       FloatValue *_distances, FloatValue *_invdistances, 
                                                       BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms, 
                              struct ResidueData *_residues, IndexValue _nres, 
                              struct ChainData *_chains, IndexValue _nchains);

#ifdef __USE_OPENCL_CODEPATH__
ErrorCode set_opencl_for_oops(cl_program *_clprogram, cl_kernel *_cldistkernel, cl_kernel *_cldistsqkernel, cl_kernel *_clinvdistkernel, cl_context *_clcontext, cl_command_queue *_clqueue, cl_mem *_clcoords, cl_mem *_cldistmask, cl_mem *_cldistances);
#endif

ErrorCode init_oops(const char* cfg);
ErrorCode finish_oops(void);
ErrorCode run_oops(void);

#endif