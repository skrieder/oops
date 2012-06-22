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
 * This file contains the definition of functions that compute torsional rotations and angles. Rotations
 * update a distance mask where each element (i, j) adopts a TRUE/FALSE value depending on whether
 * atoms with indexes i and j change their relative positions as a result of a torsional rotation. This mask
 * can be later used in the distance calculations to mask out atomic pairs.
 *
 */

#ifndef __ROTATION_H__
#define __ROTATION_H__

#include "datatypes.h"

/**
 * Puts in *angle the value of residues[res].bonds[bond].angle.
 * @param angle FloatValue*
 * @param bond IndexValue
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @return ErrorCode
 */
ErrorCode get_torsional_angle(FloatValue *angle, IndexValue bond, struct ResidueData *residues, IndexValue res);
/**
 * Sets residues[res].bonds[bond].newangle to newangle.
 * @param newangle FloatValue
 * @param bond IndexValue
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @return ErrorCode
 */
ErrorCode set_torsional_angle(FloatValue newangle, IndexValue bond, struct ResidueData *residues, IndexValue res);
/**
 * Recalculates the torsional angle of residue res using the coordinates for atom0, atom1, atom2 and atom3 given in the array coords.
 * @param bond IndexValue
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param coords FloatValue*
 * @return ErrorCode
 */
ErrorCode update_torsional_angle(IndexValue bond, struct ResidueData *residues, IndexValue res, FloatValue *coords);
/**
 * Recalculates all the torsional angles of residue res using the atomic coordinates given in the array coords.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param coords FloatValue*
 * @return ErrorCode
 */
ErrorCode update_all_torsional_angles_in_res(struct ResidueData *residues, IndexValue res, FloatValue *coords);
/**
 * Recalculates all the torsional angles in chain using the atomic coordinates given in the array coords.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param chains struct ChainData*
 * @param chain IndexValue
 * @param coords FloatValue*
 * @return ErrorCode
 */
ErrorCode update_all_torsional_angles_in_chain(struct ResidueData *residues, struct ChainData *chains, IndexValue chain, FloatValue *coords);
/**
 * Recalculates all the torsional angles in all the chains using the atomic coordinates given in the array coords.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param chains struct ChainData*
 * @param nchains IndexValue
 * @param coords FloatValue*
 * @return ErrorCode
 */
ErrorCode update_all_torsional_angles(struct ResidueData *residues, struct ChainData *chains, IndexValue nchains, FloatValue *coords);

/**
 * Updates the xyz coordinates in the array coords as the result of applying a rotation of amount rotangle (expressed in degrees) around bond in residue res.
 * @param coords FloatValue*
 * @param rotangle FloatValue
 * @param bond IndexValue
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @return ErrorCode
 */
ErrorCode rotate_coords(FloatValue *coords, FloatValue rotangle, IndexValue bond, struct ResidueData *residues, IndexValue res);

/**
 * Updates the entries of the distance mask distmask so that all the atom pairs that have their relative positions modified as result as a torsional rotation 
 * around bond in residues res are marked with TRUE in the mask.
 * @param distmask BoolValue*
 * @param atoms struct AtomData*
 * @param natoms IndexValue
 * @param bond IndexValue
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param chains struct ChainData*
 * @param nchains IndexValue
 * @return ErrorCode
 */
ErrorCode update_distmask(BoolValue *distmask, struct AtomData *atoms, IndexValue natoms, IndexValue bond, struct ResidueData *residues, IndexValue res, struct ChainData *chains, IndexValue nchains);

/**
 * Updates the xyz coordinates in coords and the distance mask distmask resulting of all the rotations ocurring at residue res. A rotation is applied when 
 * the difference between the angle and newangle values in a bond in the residue is greater than the minimum float difference allowed with the current
 * precision level. Update of the distance mask can be disabled by passing FALSE in updatemask.
 * @param res IndexValue
 * @param coords FloatValue*
 * @param distmask BoolValue*
 * @param atoms struct AtomData*
 * @param natoms IndexValue
 * @param residues struct ResidueData*
 * @param chains struct ChainData*
 * @param nchains IndexValue
 * @param createsc updatemask
 * @return ErrorCode
 */
ErrorCode update_res_coords(IndexValue res, FloatValue *coords, BoolValue *distmask, struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, struct ChainData *chains, IndexValue nchains, BoolValue updatemask);

/**
 * Updates the xyz coordinates in coords and the distance mask distmask resulting of all the rotations ocurring at chain. A rotation is applied when 
 * the difference between the angle and newangle values in a bond in the residue is greater than the minimum float difference allowed with the current
 * precision level. Update of the distance mask can be disabled by passing FALSE in updatemask.
 * @param chain IndexValue
 * @param coords FloatValue*
 * @param distmask BoolValue*
 * @param atoms struct AtomData*
 * @param natoms IndexValue
 * @param residues struct ResidueData*
 * @param chains struct ChainData*
 * @param nchains IndexValue
 * @param createsc updatemask
 * @return ErrorCode
 */
ErrorCode update_chain_coords(IndexValue chain, FloatValue *coords, BoolValue *distmask, struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, struct ChainData *chains, IndexValue nchains, BoolValue updatemask);

/**
 * Updates the xyz coordinates in coords and the distance mask distmask resulting of all the rotations ocurring across all the chains. A rotation is applied 
 * when  the difference between the angle and newangle values in a bond in the residue is greater than the minimum float difference allowed with the 
 * current precision level. Update of the distance mask can be disabled by passing FALSE in updatemask.
 * @param coords FloatValue*
 * @param distmask BoolValue*
 * @param atoms struct AtomData*
 * @param natoms IndexValue
 * @param residues struct ResidueData*
 * @param nres IndexValue
 * @param chains struct ChainData*
 * @param nchains IndexValue
 * @param createsc updatemask
 * @return ErrorCode
 */
ErrorCode update_all_coords(FloatValue *coords, BoolValue *distmask, struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, IndexValue nres, struct ChainData *chains, IndexValue nchains, BoolValue updatemask);

#endif