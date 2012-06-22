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
 * This file contains the implementation of functions that compute distance between atom pairs.
 * The distance calculations can be masked using the distance mask, so that pairs that
 * have their corresponding mask value set to FALSE don't have their distance evaluated.
 * This mask can be generated with the rotation routines defined in rotation.h.
 *
 */

#include <math.h>
#include "memory.h"
#include "distance.h"

ErrorCode clear_distmask(BoolValue *distmask, IndexValue natoms)
{
    return clear_bool_array(distmask, natoms * natoms);
}

ErrorCode calculate_all_distances(FloatValue *distances, FloatValue *coords, IndexValue natoms)
{
    IndexValue i, j;
    FloatValue xi, yi, zi;
    FloatValue xj, yj, zj;
    FloatValue dx, dy, dz;
    FloatValue dist;
    for (i = 0; i < natoms; i++)
    {
        xi = coords[xcoord_index(i)];
        yi = coords[ycoord_index(i)];
        zi = coords[zcoord_index(i)];
        for (j = i + 1; j < natoms; j++)
        {
            xj = coords[xcoord_index(j)];
            yj = coords[ycoord_index(j)];
            zj = coords[zcoord_index(j)];
            
            dx = xi - xj;
            dy = yi - yj;
            dz = zi - zj;

            #ifdef __USE_SQUARED_DISTANCES__
            dist = dx * dx + dy * dy + dz * dz;
            #else
            dist = sqrt(dx * dx + dy * dy + dz * dz);
            #endif
            
            distances[distance_index(i, j, natoms)] = dist;
            distances[distance_index(j, i, natoms)] = dist;
        }
    }
    return NO_ERROR;
}

ErrorCode calculate_masked_distances(FloatValue *distances, FloatValue *coords, IndexValue natoms, BoolValue *distmask)
{
    IndexValue i, j;
    FloatValue xi, yi, zi;
    FloatValue xj, yj, zj;
    FloatValue dx, dy, dz;
    FloatValue dist;
    for (i = 0; i < natoms; i++)
    {
        xi = coords[xcoord_index(i)];
        yi = coords[ycoord_index(i)];
        zi = coords[zcoord_index(i)];
        for (j = i + 1; j < natoms; j++)
            if (distmask[distance_index(i, j, natoms)] == 1)
            {
                xj = coords[xcoord_index(j)];
                yj = coords[ycoord_index(j)];
                zj = coords[zcoord_index(j)];

                dx = xi - xj;
                dy = yi - yj;
                dz = zi - zj;

                #ifdef __USE_SQUARED_DISTANCES__
                dist = dx * dx + dy * dy + dz * dz;
                #else
                dist = sqrt(dx * dx + dy * dy + dz * dz);
                #endif
            
                distances[distance_index(i, j, natoms)] = dist;
                distances[distance_index(j, i, natoms)] = dist;
            }
    }
    return NO_ERROR;
}

ErrorCode calculate_all_invdistances(FloatValue *invdistances, FloatValue *coords, IndexValue natoms)
{
    IndexValue i, j;
    FloatValue xi, yi, zi;
    FloatValue xj, yj, zj;
    FloatValue dx, dy, dz;
    FloatValue dist;
    for (i = 0; i < natoms; i++)
    {
        xi = coords[xcoord_index(i)];
        yi = coords[ycoord_index(i)];
        zi = coords[zcoord_index(i)];
        for (j = i + 1; j < natoms; j++)
        {
            xj = coords[xcoord_index(j)];
            yj = coords[ycoord_index(j)];
            zj = coords[zcoord_index(j)];
            
            dx = xi - xj;
            dy = yi - yj;
            dz = zi - zj;

            #ifdef __USE_SQUARED_DISTANCES__
            dist = dx * dx + dy * dy + dz * dz;
            #else
            dist = sqrt(dx * dx + dy * dy + dz * dz);
            #endif
            
            invdistances[distance_index(i, j, natoms)] = 1.0 / dist;
            invdistances[distance_index(j, i, natoms)] = 1.0 / dist;
        }
    }
    return NO_ERROR;
}

ErrorCode calculate_masked_invdistances(FloatValue *invdistances, FloatValue *coords, IndexValue natoms, BoolValue *distmask)
{
    IndexValue i, j;
    FloatValue xi, yi, zi;
    FloatValue xj, yj, zj;
    FloatValue dx, dy, dz;
    FloatValue dist;
    for (i = 0; i < natoms; i++)
    {
        xi = coords[xcoord_index(i)];
        yi = coords[ycoord_index(i)];
        zi = coords[zcoord_index(i)];
        for (j = i + 1; j < natoms; j++)
            if (distmask[distance_index(i, j, natoms)] == 1)
            {
                xj = coords[xcoord_index(j)];
                yj = coords[ycoord_index(j)];
                zj = coords[zcoord_index(j)];

                dx = xi - xj;
                dy = yi - yj;
                dz = zi - zj;

                #ifdef __USE_SQUARED_DISTANCES__
                dist = dx * dx + dy * dy + dz * dz;
                #else
                dist = sqrt(dx * dx + dy * dy + dz * dz);
                #endif
            
                invdistances[distance_index(i, j, natoms)] = 1.0 / dist;
                invdistances[distance_index(j, i, natoms)] = 1.0 / dist;
            }
    }
    return NO_ERROR;
}
