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
 * This file contains the definition of the routines to initialize and dispose the data
 *  structures.
 *
 */

#ifndef __MEMORY_H__
#define __MEMORY_H__

#include "datatypes.h"

/**
 * Allocates memory for a floating point array farray with nelements.
 * @param farray FloatValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode create_float_array(FloatValue **farray, IndexValue nelements);

/**
 * Frees the memory allocated for the floating point array farray,
 * @param farray FloatValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode delete_float_array(FloatValue **farray, IndexValue nelements);

/**
 * Allocates memory for an integer array iarray with nelements.
 * @param iarray IntValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode create_int_array(IntValue **iarray, IndexValue nelements);
/**
 * Frees the memory allocated for the integer array iarray,
 * @param iarray IntValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode delete_int_array(IntValue **iarray, IndexValue nelements);

/**
 * Allocates memory for an index array idxarray with nelements.
 * @param idxarray IndexValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode create_index_array(IndexValue **idxarray, IndexValue nelements);
/**
 * Frees the memory allocated for the index array idxarray,
 * @param idxarray IndexValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode delete_index_array(IndexValue **idxarray, IndexValue nelements);

/**
 * Allocates memory for a boolean array barray with nelements.
 * @param barray BoolValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode create_bool_array(BoolValue **barray, IndexValue nelements);
/**
 * Frees the memory allocated for the boolean array barray,
 * @param barray BoolValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode delete_bool_array(BoolValue **barray, IndexValue nelements);

/**
 * Allocates memory for a xywz coordinates array for natoms. This just means creating a floating-point
 * array with 4 * natoms.
 * @param coords FloatValue**
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode create_coords_array(FloatValue **coords, IndexValue natoms);
/**
 * Frees the memory allocated for the coordinates array coords,
 * @param coords FloatValue**
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode delete_coords_array(FloatValue **coords, IndexValue natoms);

/**
 * Allocates memory for a distance array for natoms. This just means creating a floating-point
 * array with natoms * natoms.
 * @param distances FloatValue**
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode create_distances_array(FloatValue **distances, IndexValue natoms);
/**
 * Frees the memory allocated for the distances array,
 * @param distances FloatValue**
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode delete_distances_array(FloatValue **distances, IndexValue natoms);

/**
 * Allocates memory for a distance mask array for natoms. This just means creating a boolean array with natoms * natoms.
 * @param mask BoolValue**
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode create_distmask_array(BoolValue **mask, IndexValue natoms);
/**
 * Frees the memory allocated for the distances mask array,
 * @param mask BoolValue**
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode delete_distmask_array(BoolValue **mask, IndexValue natoms);

/**
 * Sets to zero all the elements in the floating-point array farray.
 * @param farray FloatValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode clear_float_array(FloatValue *farray, IndexValue nelements);
/**
 * Sets to zero all the elements in the integer array iarray.
 * @param iarray IntValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode clear_int_array(IntValue *iarray, IndexValue nelements);
/**
 * Sets to zero all the elements in the index array idxarray.
 * @param idxarray IndexValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode clear_index_array(IndexValue *idxarray, IndexValue nelements);
/**
 * Sets to zero all the elements in the boolean array barray.
 * @param barray BoolValue**
 * @param nelements IndexValue
 * @return ErrorCode
 */
ErrorCode clear_bool_array(BoolValue *barray, IndexValue nelements);

/**
 * Allocates memory for an atoms array with natoms elements.
 * @param atoms struct AtomData**
 * @param natoms IndexValue
 * @return ErrorCode
 */
ErrorCode create_atoms_array(struct AtomData **atoms, IndexValue natoms);
/**
 * Frees the memory allocated for the atoms array, The natoms argument is set to zero.
 * @param atoms struct AtomData**
 * @param natoms IndexValue*
 * @return ErrorCode
 */
ErrorCode delete_atoms_array(struct AtomData **atoms, IndexValue *natoms);

/**
 * Allocates memory for a bonds array with nbonds elements.
 * @param bonds struct BondData**
 * @param nbonds IndexValue
 * @return ErrorCode
 */
ErrorCode create_bonds_array(struct BondData **bonds, IndexValue nbonds);
/**
 * Frees the memory allocated for the bonds array, The nbonds argument is set to zero.
 * @param bonds struct BondData**
 * @param nbonds IndexValue*
 * @return ErrorCode
 */
ErrorCode delete_bonds_array(struct BondData **bonds, IndexValue *nbonds);

/**
 * Allocates memory for a residues array with nres elements.
 * @param residues struct ResidueData**
 * @param nres IndexValue
 * @return ErrorCode
 */
ErrorCode create_residues_array(struct ResidueData **residues, IndexValue nres);
/**
 * Frees the memory allocated for the residues array, The nres argument is set to zero.
 * @param residues struct ResidueData**
 * @param nres IndexValue*
 * @return ErrorCode
 */
ErrorCode delete_residues_array(struct ResidueData **residues, IndexValue *nres);

/**
 * Allocates memory for a chains array with nchains elements.
 * @param chains struct ChainData**
 * @param nchains IndexValue
 * @return ErrorCode
 */
ErrorCode create_chains_array(struct ChainData **chains, IndexValue nchains);
/**
 * Frees the memory allocated for the chains array, The nchains argument is set to zero.
 * @param chains struct ChainData**
 * @param nchains IndexValue*
 * @return ErrorCode
 */
ErrorCode delete_chains_array(struct ChainData **chains, IndexValue *nchains);

#endif