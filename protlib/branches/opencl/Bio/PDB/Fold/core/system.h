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
 * Global system variables that store molecule data. These variables are shared accross the entire library.
 *
 */

#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#include "datatypes.h"

#ifdef __USE_OPENCL_CODEPATH__
#include <CL/cl.h>

// OpenCL required variables.
static cl_context clcontext = 0;
static cl_device_id *cldevices = 0;

static cl_command_queue clqueue = 0;
static cl_mem clcoords = 0;
static cl_mem cldistmask = 0;
static cl_mem cldistances = 0;

static cl_program clprogram = 0;
static cl_kernel cldistkernel = 0;
static cl_kernel cldistsqkernel = 0;
static cl_kernel clinvdistkernel = 0;
#endif

static struct AtomData *atoms;          // Array of atom data.
static FloatValue *coords;              // Array of atom coordinates.
static FloatValue *distances;           // Array of atom pair-wise distances.
static FloatValue *invdistances;        // Array of atom pair-wise 1/distances.
static BoolValue *atommask;             // Array of atom masking values.
static BoolValue *distmask;             // Array of distance masking values.  
static IndexValue natoms;               // Total number of atoms.

static struct ResidueData *residues;    // Array of residue data.
static IndexValue nres;                 // Total number of residues.
 
static struct ChainData *chains;        // Array of chain data.
static IndexValue nchains;              // Total number of chains.

// Convenience define to access the elements of the distance matrix (which is stored as a 1D array).
#define distance_index_natoms(atom0, atom1)     (atom0 * natoms + atom1)

#endif
