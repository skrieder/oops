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
 * This file contains the implementation of functions that compute torsional rotations and angles. Rotations
 * update a distance mask where each element (i, j) adopts a TRUE/FALSE value depending on whether
 * atoms with indexes i and j change their relative positions as a result of a torsional rotation. This mask
 * can be later used in the distance calculations to mask out atomic pairs.
 *
 */

#include <math.h>
#include "vector.h"
#include "datautils.h"
#include "rotation.h"

ErrorCode get_torsional_angle(FloatValue *angle, IndexValue bond, struct ResidueData *residues, IndexValue res)
{
    if (residues[res].bonds[bond].type == NOT_SET)
    {
        return ACCESSING_UNSET_BOND_ERROR;
    }

    *angle = residues[res].bonds[bond].angle;

    return NO_ERROR;
}

ErrorCode set_torsional_angle(FloatValue newangle, IndexValue bond, struct ResidueData *residues, IndexValue res)
{
    if (residues[res].bonds[bond].type == NOT_SET)
    {
        return ACCESSING_UNSET_BOND_ERROR;
    }

    residues[res].bonds[bond].newangle = newangle;

    return NO_ERROR;
}

ErrorCode update_torsional_angle(IndexValue bond, struct ResidueData *residues, IndexValue res, FloatValue *coords)
{
    if (residues[res].bonds[bond].type == NOT_SET)
    {
        return ACCESSING_UNSET_BOND_ERROR;
    }

    IndexValue atom0, atom1, atom2, atom3;
    FloatValue x0, y0, z0;
    FloatValue x1, y1, z1;
    FloatValue x2, y2, z2;
    FloatValue x3, y3, z3;

    FloatValue x01, y01, z01;
    FloatValue x12, y12, z12;
    FloatValue x32, y32, z32;

    FloatValue px, py, pz;
    FloatValue qx, qy, qz;
    FloatValue rx, ry, rz;

    FloatValue u, v, u1, v1;
    FloatValue a;

    atom0 = residues[res].bonds[bond].atom0;
    x0 = coords[xcoord_index(atom0)];
    y0 = coords[ycoord_index(atom0)];
    z0 = coords[zcoord_index(atom0)];

    atom1 = residues[res].bonds[bond].atom1;
    x1 = coords[xcoord_index(atom1)];
    y1 = coords[ycoord_index(atom1)];
    z1 = coords[zcoord_index(atom1)];
    
    atom2 = residues[res].bonds[bond].atom2;
    x2 = coords[xcoord_index(atom2)];
    y2 = coords[ycoord_index(atom2)];
    z2 = coords[zcoord_index(atom2)];

    atom3 = residues[res].bonds[bond].atom3;
    x3 = coords[xcoord_index(atom3)];
    y3 = coords[ycoord_index(atom3)];
    z3 = coords[zcoord_index(atom3)];

    vector_subtract(&x01, &y01, &z01, x0, y0, z0, x1, y1, z1);     // v01 = v0 - v1
    vector_subtract(&x32, &y32, &z32, x3, y3, z3, x2, y2, z2);     // v32 = v3 - v2
    vector_subtract(&x12, &y12, &z12, x1, y1, z1, x2, y2, z2);     // v12 = v1 - v2

    vector_crossprod(&px, &py, &pz, x12, y12, z12, x01, y01, z01); // p = v12 x v01
    vector_crossprod(&qx, &qy, &qz, x12, y12, z12, x32, y32, z32); // q = v12 x v32
    vector_crossprod(&rx, &ry, &rz, x12, y12, z12, qx, qy, qz);    // r = v12 x q

    vector_dotprod(&u, qx, qy, qz, qx, qy, qz);                    // u = q * q
    vector_dotprod(&v, rx, ry, rz, rx, ry, rz);                    // v = r * r

    if (u <= 0.0 || v <= 0.0)
    {
        a = 360.0;
    }
    else 
    {
        vector_dotprod(&u1, px, py, pz, qx, qy, qz);               // u1 = p * q
        vector_dotprod(&v1, px, py, pz, rx, ry, rz);               // v1 = p * r

        u = u1 / sqrt(u);
        v = v1 / sqrt(v);

        if (fabs(u) > MIN_FLOAT_DIFF || fabs(v) > MIN_FLOAT_DIFF) a = atan2(v, u) * RADIANS_TO_DEGREES;
        else a = 360.0;
    }

    residues[res].bonds[bond].angle = a;
    residues[res].bonds[bond].newangle = a;

    return NO_ERROR;
}

ErrorCode update_all_torsional_angles_in_res(struct ResidueData *residues, IndexValue res, FloatValue *coords)
{
    IndexValue bond;
    ErrorCode error0, error;
    error0 = NO_ERROR;
    for (bond = 0; bond < residues[bond].nbonds; bond++)
    {
        error = update_torsional_angle(bond, residues, res, coords);
        if (error != NO_ERROR) error0 = error;
    }
    return error0;
}

ErrorCode update_all_torsional_angles_in_chain(struct ResidueData *residues, struct ChainData *chains, IndexValue chain, FloatValue *coords)
{
    IndexValue res;
    ErrorCode error0, error;
    error0 = NO_ERROR;
    for (res = chains[chain].firstres; res <= chains[chain].lastres; res++)
    {
        error = update_all_torsional_angles_in_res(residues, res, coords);
        if (error != NO_ERROR) error0 = error;
    }
    return error0;
}

ErrorCode update_all_torsional_angles(struct ResidueData *residues, struct ChainData *chains, IndexValue nchains, FloatValue *coords)
{
    IndexValue chain;
    ErrorCode error0, error;
    error0 = NO_ERROR;
    for (chain = 0; chain < nchains; chain++)
    {
        error = update_all_torsional_angles_in_chain(residues, chains, chain, coords);
        if (error != NO_ERROR) error0 = error;
    }
    return error0;
}

ErrorCode rotate_coords(FloatValue *coords, FloatValue rotangle, IndexValue bond, struct ResidueData *residues, IndexValue res)
{
    IndexValue atom1, atom2;
    FloatValue x1, y1, z1;
    FloatValue x2, y2, z2;

    FloatValue rx, ry, rz;

    FloatValue matxx, matxy, matxz;
    FloatValue matyx, matyy, matyz;
    FloatValue matzx, matzy, matzz;

    IndexValue lastatom;

    FloatValue dx, dy, dz;
    FloatValue x, y, z;
    IndexValue atom;

    atom1 = residues[res].bonds[bond].atom1;
    x1 = coords[xcoord_index(atom1)];
    y1 = coords[ycoord_index(atom1)];
    z1 = coords[zcoord_index(atom1)];
    
    atom2 = residues[res].bonds[bond].atom2;
    x2 = coords[xcoord_index(atom2)];
    y2 = coords[ycoord_index(atom2)];
    z2 = coords[zcoord_index(atom2)];
 
    // This is the index of the last atom in the array affected by this rotation,.
    lastatom = residues[res].bonds[bond].lastatom;

    // Calculating normalized rotation axis.
    vector_subtract(&rx, &ry, &rz, x2, y2, z2, x1, y1, z1);
    vector_normalize(&rx, &ry, &rz);

    // Calculating rotation matrix for rotation angle rotangle around rotation axis (rx, ry, rz).
    create_rotation_matrix(&matxx, &matxy, &matxz, &matyx, &matyy, &matyz, &matzx, &matzy, &matzz, rotangle, rx, ry, rz);

    for (atom = atom2 + 1; atom <= lastatom; atom++) //previous version had atom<lastatom. stranded last atom in space
    {
        dx = coords[xcoord_index(atom)] - x2;
        dy = coords[ycoord_index(atom)] - y2;
        dz = coords[zcoord_index(atom)] - z2;

        x = matxx * dx + matxy * dy + matxz * dz;
        y = matyx * dx + matyy * dy + matyz * dz;
        z = matzx * dx + matzy * dy + matzz * dz;

        coords[xcoord_index(atom)] = x + x2; 
        coords[ycoord_index(atom)] = y + y2; 
        coords[zcoord_index(atom)] = z + z2;
    }

    return NO_ERROR;
}

ErrorCode update_distmask(BoolValue *distmask, struct AtomData *atoms, IndexValue natoms, IndexValue bond, struct ResidueData *residues, IndexValue res, struct ChainData *chains, IndexValue nchains)
{
    IndexValue atom1, atom2;
    IndexValue bondlastatom;

    IndexValue chain, c;
    IndexValue firstres, lastres;
    IndexValue firstatom, lastatom;

    // These two atoms define the rotating bond.
    atom1 = residues[res].bonds[bond].atom1;
    atom2 = residues[res].bonds[bond].atom2;

    // And this is the last atom in the chain affected by the rotation at bond.
    bondlastatom = residues[res].bonds[bond].lastatom;

    // Getting first and last atom in chain.
    chain = residues[res].chain;
    firstres = chains[chain].firstres;
    lastres = chains[chain].lastres;
    firstatom = residues[firstres].firstatom;
    lastatom = residues[lastres].lastatom;

    // The atoms in the interval [atom2, bondlastatom], affected by the rotation at bond in res, change their relative positions to the atoms within
    // all the previous atoms in the chain (intervals [firstatom, atom1],and [atom1 + 1, atom2 - 1]) and all the subsequent atoms in the chain, which 
    // are contained in the interval [atom2 + 1, lastatom].
    set_distmask_intervals_product(1, distmask, atoms, natoms, firstatom, atom1, atom2, bondlastatom);
    set_distmask_intervals_product(1, distmask, atoms, natoms, atom1 + 1, atom2 - 1, atom2, bondlastatom);
    set_distmask_intervals_product(1, distmask, atoms, natoms, atom2 + 1, lastatom, atom2, bondlastatom);

    // The rotated interval [atom2, bondlastatom], has also changed its position with respect to all the other chains.
    for (c = 0; c < nchains; c++)
        if (c != chain)
        {
            // Getting first and last atom in chain c.
            firstres = chains[c].firstres;
            lastres = chains[c].lastres;
            firstatom = residues[firstres].firstatom;
            lastatom = residues[lastres].lastatom;
            set_distmask_intervals_product(1, distmask, atoms, natoms, firstatom, lastatom, atom2, bondlastatom);
        }

    return NO_ERROR;
}

ErrorCode update_res_coords(IndexValue res, FloatValue *coords, BoolValue *distmask, struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, struct ChainData *chains, IndexValue nchains, BoolValue updatemask)
{
    IndexValue bond;
    FloatValue alpha;
    for (bond = 0; bond < residues[res].nbonds; bond++)
        if (residues[res].bonds[bond].type != NOT_SET)
        {
            alpha = residues[res].bonds[bond].newangle - residues[res].bonds[bond].angle;
	    //printf("res: %i alpha: %f newangle %f oldangle %f\n",res,alpha, residues[res].bonds[bond].newangle,residues[res].bonds[bond].angle);
            if (fabs(alpha) > MIN_FLOAT_DIFF)
            {
	       
                //rotate_coords(coords, alpha, bond, residues, res);
                rotate_coords(coords, -alpha, bond, residues, res);
                
                if (updatemask)
                {
                    update_distmask(distmask, atoms, natoms, bond, residues, res, chains, nchains);
                }


                residues[res].bonds[bond].angle += alpha; //Aashish: change this to residues[res].bonds[bond].angle=residues[res].bonds[bond].newangle ?

	    //printf(" Getting res: %i alpha: %f newangle %f oldangle %f\n",res,alpha, residues[res].bonds[bond].newangle,residues[res].bonds[bond].angle);

	    }
        }
    return NO_ERROR;
}

ErrorCode update_chain_coords(IndexValue chain, FloatValue *coords, BoolValue *distmask, struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, struct ChainData *chains, IndexValue nchains, BoolValue updatemask)
{
    IndexValue firstres = chains[chain].firstres;
    IndexValue lastres = chains[chain].lastres;
    IndexValue res;
    for (res = firstres; res <= lastres; res++)
    {
        update_res_coords(res, coords, distmask, atoms, natoms, residues, chains, nchains, updatemask);
    }
    return NO_ERROR;
}

ErrorCode update_all_coords(FloatValue *coords, BoolValue *distmask, struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, IndexValue nres, struct ChainData *chains, IndexValue nchains, BoolValue updatemask)
{
    IndexValue res;
    for (res = 0; res < nres; res++)
    {
        update_res_coords(res, coords, distmask, atoms, natoms, residues, chains, nchains, updatemask);
    }
    return NO_ERROR;
}
