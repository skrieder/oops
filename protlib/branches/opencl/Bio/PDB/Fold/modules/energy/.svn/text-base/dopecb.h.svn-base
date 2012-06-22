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
 % WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

/**
 * DOPE-CB statisical potential function.
 * Authors: Min-yi Shen (original creator of DOPE), Andres Colubri (protlib implementation), Glen Hocky (protlib2 port)
 *
 */

#ifndef __DOPE_H__
#define __DOPE_H__

#include "../../core/datatypes.h"

ErrorCode set_system_for_dopecb(struct AtomData *_atoms, FloatValue *_coords, 
                                                         FloatValue *_distances, FloatValue *_invdistances, 
                                                         BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms, 
                                struct ResidueData *_residues, IndexValue _nres, 
                                struct ChainData *_chains, IndexValue _nchains);

ErrorCode load_dopecb_parfile(const char* parfilename);

ErrorCode init_dopecb(const char* cfg);
ErrorCode finish_dopecb(void);

ErrorCode pre_dopecb_calc(void);
ErrorCode post_dopecb_calc(void);

ErrorCode calc_dopecb_energy(FloatValue *energy);
ErrorCode calc_dopecb_gradient(FloatValue *gradient);

#define convert_to_dopecb_row(aa0,at0,aa1,at1) (aa0*38400+aa1*1920+at0*240+at1*30)
#define convert_to_dopecb_pos(aa0,at0,aa1,at1,bin) (aa0*38400+aa1*1920+at0*240+at1*30+bin)

//IntValue convert_to_dope_row(TypeCode aa0,TypeCode at0,TypeCode aa1,TypeCode at1){return aa0*38400+aa1*1920+at0*240+at1*30;}
//IntValue convert_to_dope_pos(TypeCode aa0,TypeCode at0,TypeCode aa1,TypeCode at1,IntValue bin){return aa0*38400+aa1*1920+at0*240+at1*30+bin;}

#endif
