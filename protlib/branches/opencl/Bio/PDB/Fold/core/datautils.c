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
 * This file implements a number of utility functions to handle the basic data structures..
 *
 */

#include "datautils.h"

ErrorCode set_atom_coords(struct AtomData *atoms, FloatValue *coords, IndexValue natoms)
{
    IndexValue i;
    for (i = 0; i < natoms; i++)
    {
        // Some pointer arithmetics here: atoms[i].x contains the address of the n-th element in the coords array, with n being
        // the position corresponding to the x coordinate of the i-th atom, i.e.: the address of element i * 4 + 0 (+1, 2, 3 in the 
        // case of y, z and w).         
        atoms[i].x = coords + (i * 4 + 0);
        atoms[i].y = coords + (i * 4 + 1);
        atoms[i].z = coords + (i * 4 + 2);
        atoms[i].w = coords + (i * 4 + 3);
    }

    return NO_ERROR;
}

ErrorCode get_atom_coords(FloatValue *x, FloatValue *y, FloatValue *z, struct AtomData *atoms, IndexValue atom)
{
    *x = *(atoms[atom].x);
    *y = *(atoms[atom].y);
    *z = *(atoms[atom].z);
    return NO_ERROR;
}

ErrorCode set_distmask_intervals_product(BoolValue value, BoolValue *mask, struct AtomData *atoms, IndexValue natoms, IndexValue atom00, IndexValue atom01, IndexValue atom10, IndexValue atom11)
{ 
    IndexValue i, j;
    for (i = atom00; i <= atom01; i++)
        for (j = atom10; j <= atom11; j++)
        {
            #ifdef __MASK_INTER_RESIDUE_PAIRS__
            if (atoms[i].res != atoms[j].res) mask[distance_index(i, j, natoms)] = value;
            #else
            mask[distance_index(i, j, natoms)] = value;
            #endif
        }
    return NO_ERROR;
}
