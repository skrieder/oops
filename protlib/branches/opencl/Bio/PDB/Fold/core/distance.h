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
 * This file contains the definition of functions that compute distance between atom pairs.
 * The distance calculations can be masked using the distance mask, so that pairs that
 * have their corresponding mask value set to FALSE don't have their distance evaluated.
 * This mask can be generated with the rotation routines defined in rotation.h.
 *
 */

#ifndef __DISTANCE_H__
#define __DISTANCE_H__

#include "datatypes.h"

/**
 * Sets all the entries in distmask to FALSE.
 * @param distmask BoolValue*
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode clear_distmask(BoolValue *distmask, IndexValue natoms);

/**
 * Calculates the distances between all the natoms, whose xyz coordinates are stored in the array coords.
 * @param distances FloatValue*
 *@param coords FloatValue*
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode calculate_all_distances(FloatValue *distances, FloatValue *coords, IndexValue natoms);

/**
 * Calculates the distances between the atom pairs that have TRUE in the correspoing entry in the distance mask distmask. 
 * @param distances FloatValue*
 *@param coords FloatValue*
 * @param natoms IndexValue
 *@param distmask BoolValue*
 * @return ErrorCode
 */
ErrorCode calculate_masked_distances(FloatValue *distances, FloatValue *coords, IndexValue natoms, BoolValue *distmask);

/**
 * Calculates the inverse distances between all the natoms, whose xyz coordinates are stored in the array coords.
 * @param distances FloatValue*
 *@param coords FloatValue*
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode calculate_all_invdistances(FloatValue *invdistances, FloatValue *coords, IndexValue natoms);

/**
 * Calculates the inverse distances between the atom pairs that have TRUE in the correspoing entry in the distance mask distmask. 
 * @param distances FloatValue*
 *@param coords FloatValue*
 * @param natoms IndexValue
 *@param distmask BoolValue*
 * @return ErrorCode
 */
ErrorCode calculate_masked_invdistances(FloatValue *invdistances, FloatValue *coords, IndexValue natoms, BoolValue *distmask);

#endif