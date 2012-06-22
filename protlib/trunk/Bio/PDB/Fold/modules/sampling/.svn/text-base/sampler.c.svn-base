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
 * Dummy Conformational Sampler provided as example.
 *
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../../core/memory.h"
#include "../../core/modules.h"
#include "../../core/fileutils.h"
#include "../../core/errorhandler.h"
#include "../../core/rotation.h"
#include "../../core/distance.h"
#include "../../core/random.h"
#include "../../core/strutils.h"
#include "sampler.h"

// Shared variables.
#include "../shared/variables.h"
extern BoolValue *flags;         
extern FloatValue T;

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

// Local copy of system variables.
struct AtomData *atoms;          // Array of atom data.
FloatValue *coords;              // Array of atom coordinates.
FloatValue *distances;           // Array of atom pair-wise distances.
FloatValue *invdistances;        // Array of atom pair-wise 1/distances.
BoolValue *atommask;             // Array of atom masking values.
BoolValue *distmask;             // Array of distance masking values.  
IndexValue natoms;               // Total number of atoms.
struct ResidueData *residues;    // Array of residue data.
IndexValue nres;                 // Total number of residues. 
struct ChainData *chains;        // Array of chain data.
IndexValue nchains;              // Total number of chains.

ErrorCode set_system_for_sampler(struct AtomData *_atoms, FloatValue *_coords, 
                                                      FloatValue *_distances, FloatValue *_invdistances, 
                                                      BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms, 
                                 struct ResidueData *_residues, IndexValue _nres, 
                                 struct ChainData *_chains, IndexValue _nchains)
{
    atoms = _atoms;
    coords = _coords;
    distances = _distances;
    invdistances = _invdistances;
    atommask = _atommask;
    distmask = _distmask;
    natoms = _natoms;
    residues = _residues;
    nres = _nres;
    chains = _chains;
    nchains = _nchains;

    return NO_ERROR;
}

ErrorCode init_sampler(const char* cfg)
{
    return NO_ERROR;
}

ErrorCode finish_sampler(void)
{
    return NO_ERROR;
}

ErrorCode pre_sampler_step(void)
{
    return NO_ERROR;
}

ErrorCode post_sampler_step(void)
{
    return NO_ERROR;
}

ErrorCode run_sampler_step(void)
{
    return NO_ERROR;
}
