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
 * This file defines a number of utility functions to handle the basic data structures..
 *
 */

#ifndef __DATAUTILS_H__
#define __DATAUTILS_H__

#include "datatypes.h"

/**
 * Sets the xyzw pointers in the atom records in the atoms array to the corresponding positions in the coordinates array coords.
 * @param atoms struct AtomData*
 * @param coords FloatValue*
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode set_atom_coords(struct AtomData *atoms, FloatValue *coords, IndexValue natoms);

/**
 * Puts the xyz coordinates of atoms[atom] into the (*x, *y, *z) floating-point variables passed as argument.
 * @param x FloatValue*
 * @param y FloatValue*
 * @param z FloatValue*
 * @param atoms struct AtomData*
 * @param atom IndexValue
 * @return ErrorCode
 */
ErrorCode get_atom_coords(FloatValue *x, FloatValue *y, FloatValue *z, struct AtomData *atoms, IndexValue atom);

/**
 * Sets the pairs (i, j) in the boolean distance mask to value, where (i, j) are all the pairs in the cartesian product between the intervals [atom00, atom01] and [atom10, atom11].
 * @param value BoolValue
 * @param mask BoolValue*
 * @param atoms struct AtomData*
 * @param natom IndexValue
 * @param atom00 IndexValue
 * @param atom01 IndexValue
 * @param atom10 IndexValue
 * @param atom11 IndexValue
 * @return ErrorCode
 */
ErrorCode set_distmask_intervals_product(BoolValue value, BoolValue *mask, struct AtomData *atoms, IndexValue natoms, IndexValue atom00, IndexValue atom01, IndexValue atom10, IndexValue atom11);

#endif