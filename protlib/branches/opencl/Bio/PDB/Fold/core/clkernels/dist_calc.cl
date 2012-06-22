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
 * Compute kernels for distance calculations.
 *
 */

__kernel void distsq_calc(__global float4 *coords, __global char *distmask, __global float *distances, int natoms)
{
    unsigned int i = get_global_id(0);
    unsigned int j = get_global_id(1);

    float4 vi, vj;

    vi = coords[i];
    vj = coords[j];

    int idx = i * natoms + j;

    distances[idx] = distmask[idx] * dot(vi, vj);
}

__kernel void invdist_calc(__global float4 *coords, __global char *distmask, __global float *distances, int natoms)
{
    unsigned int i = get_global_id(0);
    unsigned int j = get_global_id(1);

    float4 vi, vj;

    vi = coords[i];
    vj = coords[j];

    int idx = i * natoms + j;

    distances[idx] = distmask[idx] * 1.0 / distance(vi, vj);
}

__kernel void dist_calc(__global float4 *coords, __global char *distmask, __global float *distances, int natoms)
{
    unsigned int i = get_global_id(0);
    unsigned int j = get_global_id(1);

    float4 vi, vj;

    vi = coords[i];
    vj = coords[j];

    int idx = i * natoms + j;

    distances[idx] = distmask[idx] * distance(vi, vj);
}
