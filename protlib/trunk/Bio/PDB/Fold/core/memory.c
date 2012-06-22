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
 * This file contains the implementation of the routines to initialize and dispose the data
 *  structures.
 *
 */

#include <string.h>
#include "memory.h"

ErrorCode create_float_array(FloatValue **farray, IndexValue nelements)
{
    FloatValue *ptr = malloc(nelements * sizeof (FloatValue));

    if (ptr == NULL) 
    {
        *farray = NULL;
        return MEMORY_ALLOCATION_ERROR;
    }

    *farray = ptr;

    // Intializing array.
    memset(ptr, 0, nelements * sizeof(FloatValue));

    return NO_ERROR;
}

ErrorCode delete_float_array(FloatValue **farray, IndexValue nelements)
{
    FloatValue *ptr = *farray;

    if (ptr == NULL && 0 < nelements)
    {
        return ARRAY_INCONSISTENCY_ERROR;
    }

    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
    *farray = NULL;

    return NO_ERROR;
}

ErrorCode create_int_array(IntValue **iarray, IndexValue nelements)
{
    IntValue *ptr = malloc(nelements * sizeof (IntValue));

    if (ptr == NULL) 
    {
        *iarray = NULL;
        return MEMORY_ALLOCATION_ERROR;
    }

    *iarray = ptr;

    // Intializing array.
    memset(ptr, 0, nelements * sizeof(IntValue));

    return NO_ERROR;
}

ErrorCode delete_int_array(IntValue **iarray, IndexValue nelements)
{
    IntValue *ptr = *iarray;

    if (ptr == NULL && 0 < nelements)
    {
        return ARRAY_INCONSISTENCY_ERROR;
    }

    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
    *iarray = NULL;

    return NO_ERROR;
}

ErrorCode create_index_array(IndexValue **idxarray, IndexValue nelements)
{
    IndexValue *ptr = malloc(nelements * sizeof (IndexValue));

    if (ptr == NULL) 
    {
        *idxarray = NULL;
        return MEMORY_ALLOCATION_ERROR;
    }

    *idxarray = ptr;

    // Intializing array.
    memset(ptr, 0, nelements * sizeof(IndexValue));

    return NO_ERROR;
}

ErrorCode delete_index_array(IndexValue **idxarray, IndexValue nelements)
{
    IndexValue *ptr = *idxarray;

    if (ptr == NULL && 0 < nelements)
    {
        return ARRAY_INCONSISTENCY_ERROR;
    }

    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
    *idxarray = NULL;

    return NO_ERROR;
}

ErrorCode create_bool_array(BoolValue **barray, IndexValue nelements)
{
    BoolValue *ptr = malloc(nelements * sizeof (BoolValue));

    if (ptr == NULL) 
    {
        *barray = NULL;
        return MEMORY_ALLOCATION_ERROR;
    }

    *barray = ptr;

    // Intializing array.
    memset(ptr, 0, nelements * sizeof(BoolValue));

    return NO_ERROR;
}

ErrorCode delete_bool_array(BoolValue **barray, IndexValue nelements)
{
    BoolValue *ptr = *barray;

    if (ptr == NULL && 0 < nelements)
    {
        return ARRAY_INCONSISTENCY_ERROR;
    }

    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
    *barray = NULL;

    return NO_ERROR;
}

ErrorCode create_coords_array(FloatValue **coords, IndexValue natoms)
{
    return create_float_array(coords, 4 * natoms);
}

ErrorCode delete_coords_array(FloatValue **coords, IndexValue natoms)
{
    return delete_float_array(coords, 4 * natoms);
}

ErrorCode create_distances_array(FloatValue **distances, IndexValue natoms)
{
    return create_float_array(distances, natoms * natoms);
}

ErrorCode delete_distances_array(FloatValue **distances, IndexValue natoms)
{
    return delete_float_array(distances, natoms * natoms);
}

ErrorCode create_distmask_array(BoolValue **distmask, IndexValue natoms)
{
    return create_bool_array(distmask, natoms * natoms);
}

ErrorCode delete_distmask_array(BoolValue **distmask, IndexValue natoms)
{
    return delete_bool_array(distmask, natoms * natoms);
}

ErrorCode clear_float_array(FloatValue *farray, IndexValue nelements)
{
    memset(farray, 0, nelements * sizeof(FloatValue));
    return NO_ERROR;
}

ErrorCode clear_int_array(IntValue *iarray, IndexValue nelements)
{
    memset(iarray, 0, nelements * sizeof(IntValue));
    return NO_ERROR;
}

ErrorCode clear_index_array(IndexValue *idxarray, IndexValue nelements)
{
    memset(idxarray, 0, nelements * sizeof(IndexValue));
    return NO_ERROR;
}

ErrorCode clear_bool_array(BoolValue *barray, IndexValue nelements)
{
    memset(barray, 0, nelements * sizeof(BoolValue));
    return NO_ERROR;
}

ErrorCode create_atoms_array(struct AtomData **atoms, IndexValue natoms)
{
    struct AtomData *ptr = malloc(natoms * sizeof (struct AtomData));

    if (ptr == NULL) 
    {
        *atoms = NULL;
        return MEMORY_ALLOCATION_ERROR;
    }

    *atoms = ptr;

    // Intializing array.
    IndexValue i;
    for (i = 0; i < natoms; i++)
    {
        ptr[i].type = NOT_SET;

        ptr[i].x = NULL;
        ptr[i].y = NULL;
        ptr[i].z = NULL;
        ptr[i].w = NULL;

        ptr[i].res = 0; 
        ptr[i].chain = 0; 

        ptr[i].niprop = 0;
        ptr[i].iproperties = NULL;
        ptr[i].nfprop = 0;
        ptr[i].fproperties = NULL;
    }

    return NO_ERROR;
}

ErrorCode delete_atoms_array(struct AtomData **atoms, IndexValue *natoms)
{
    struct AtomData *ptr = *atoms;

    if (ptr == NULL && 0 < *natoms)
    {
        return ARRAY_INCONSISTENCY_ERROR;
    }

    // Releasing memory of arrays (bonds, iproperties, fproperties) inside the residues array.
    IndexValue i;
    for (i = 0; i < *natoms; i++)
    {
        if (ptr[i].iproperties != NULL) free(ptr[i].iproperties);
        if (ptr[i].fproperties != NULL) free(ptr[i].fproperties);
    }

    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
    *atoms = NULL;
    *natoms = 0;

    return NO_ERROR;
}

ErrorCode create_bonds_array(struct BondData **bonds, IndexValue nbonds)
{
    struct BondData *ptr = malloc(nbonds * sizeof (struct BondData));

    if (ptr == NULL) 
    {
        *bonds = NULL;
        return MEMORY_ALLOCATION_ERROR;
    }

    *bonds = ptr;
    
    // Intializing array.
    IndexValue i;
    for (i = 0; i < nbonds; i++)
    {
        ptr[i].type = NOT_SET;

        ptr[i].atom0 = 0;
        ptr[i].atom1 = 0;
        ptr[i].atom2 = 0;
        ptr[i].atom3 = 0;

        ptr[i].lastatom = 0;

        ptr[i].angle = 0.0;
        ptr[i].newangle = 0.0;
    }

    return NO_ERROR;
}

ErrorCode delete_bonds_array(struct BondData **bonds, IndexValue *nbonds)
{
    struct BondData *ptr = *bonds;

    if (ptr == NULL && 0 < *nbonds)
    {
        return ARRAY_INCONSISTENCY_ERROR;
    }

    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
    *bonds = NULL;
    *nbonds = 0;

    return NO_ERROR;
}

ErrorCode create_residues_array(struct ResidueData **residues, IndexValue nres)
{
    struct ResidueData *ptr = malloc(nres * sizeof (struct ResidueData));

    if (ptr == NULL) 
    {
        *residues = NULL;
        return MEMORY_ALLOCATION_ERROR;
    }

    *residues = ptr;

    // Intializing array.
    IndexValue i;
    for (i = 0; i < nres; i++)
    {
        ptr[i].type = NOT_SET;

        ptr[i].firstatom = 0;
        ptr[i].lastatom = 0;

        ptr[i].chain = 0; 

        ptr[i].nbonds = 0;
        ptr[i].bonds = NULL;

        ptr[i].niprop = 0;
        ptr[i].iproperties = NULL;
        ptr[i].nfprop = 0;
        ptr[i].fproperties = NULL;
    }

    return NO_ERROR;
}

ErrorCode delete_residues_array(struct ResidueData **residues, IndexValue *nres)
{
    struct ResidueData *ptr = *residues;

    if (ptr == NULL && 0 < *nres)
    {
        return ARRAY_INCONSISTENCY_ERROR;
    }

    // Releasing memory of arrays (bonds, iproperties, fproperties) inside the residues array.
    IndexValue i;
    for (i = 0; i < *nres; i++)
    {
        if (ptr[i].bonds != NULL) free(ptr[i].bonds);
        if (ptr[i].iproperties != NULL) free(ptr[i].iproperties);
        if (ptr[i].fproperties != NULL) free(ptr[i].fproperties);
    }

    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
    *residues = NULL;
    *nres = 0;

    return NO_ERROR;
}

ErrorCode create_chains_array(struct ChainData **chains, IndexValue nchains)
{
    struct ChainData *ptr = malloc(nchains * sizeof (struct ChainData));

    if (ptr == NULL) 
    {
        *chains = NULL;
        return MEMORY_ALLOCATION_ERROR;
    }

    *chains = ptr;

    // Intializing array.
    IndexValue i;
    for (i = 0; i < nchains; i++)
    {
        ptr[i].type = NOT_SET;

        ptr[i].firstres = 0;
        ptr[i].lastres = 0;
    }

    return NO_ERROR;
}

ErrorCode delete_chains_array(struct ChainData **chains, IndexValue *nchains)
{
    struct ChainData *ptr = *chains;

    if (ptr == NULL && 0 < *nchains)
    {
        return ARRAY_INCONSISTENCY_ERROR;
    }

    if (ptr != NULL)
    {
        free(ptr);
        ptr = NULL;
    }
    *chains = NULL;
    *nchains = 0;

    return NO_ERROR;
}
