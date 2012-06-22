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
 * Dummy Energy Function provided as example.
 *
 */

#ifndef __ENERGY_H__
#define __ENERGY_H__

#include "../../core/datatypes.h"

ErrorCode set_system_for_energy(struct AtomData *_atoms, FloatValue *_coords, 
                                                       FloatValue *_distances, FloatValue *_invdistances, 
                                                       BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms, 
                                struct ResidueData *_residues, IndexValue _nres, 
                                struct ChainData *_chains, IndexValue _nchains);

ErrorCode init_energy(const char* cfg);
ErrorCode finish_energy(void);

ErrorCode pre_energy_calc(void);
ErrorCode post_energy_calc(void);

ErrorCode calc_energy_energy(FloatValue *energy);
ErrorCode calc_energy_gradient(FloatValue *gradient);

#endif
