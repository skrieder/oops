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
 * This file implements various debugging and testing routines. New user-defined debug functions
 *. should be added here.
 *
 */

#include <stdio.h>
#include "debugutils.h"

ErrorCode print_molecule_data(struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, IndexValue nres, struct ChainData *chains, IndexValue nchains)
{
    IndexValue i, j;
    TypeCode type;
    IndexValue atom0, atom1, atom2, atom3, lastatom;
    FloatValue angle;

    printf("Atom data:\n");
    for (i = 0; i < natoms; i++)
    {
        printf("Atom %d type:%d coords:[%f,%f,%f] res:%d chain:%d\n", i, atoms[i].type, *(atoms[i].x), *(atoms[i].y), *(atoms[i].z), atoms[i].res, atoms[i].chain);
    }

    printf("Residue data:\n");
    for (i = 0; i < nres; i++)
    {
        printf("Residue %d type:%d firstatom:%d lastatom:%d chain:%d nbonds:%d\n", i, residues[i].type, residues[i].firstatom, residues[i].lastatom, atoms[i].chain, residues[i].nbonds);
        printf("    Bond data:\n"); 
        for (j = 0; j < residues[i].nbonds; j++)
        {
            type = residues[i].bonds[j].type;
            atom0 = residues[i].bonds[j].atom0;
            atom1 = residues[i].bonds[j].atom1;
            atom2 = residues[i].bonds[j].atom2;
            atom3 = residues[i].bonds[j].atom3;
            lastatom = residues[i].bonds[j].lastatom;
            angle = residues[i].bonds[j].angle;
            printf("    bond %d type:%d atom0:%d[type:%d] atom1:%d[type:%d]\n", i, type, atom0, atoms[atom0].type, atom1, atoms[atom1].type);
            printf("                    atom2:%d[type:%d] atom3:%d[type:%d]\n", atom2, atoms[atom2].type, atom3, atoms[atom3].type);
            printf("                    lastatom:%d[type:%d] angle:%f\n", lastatom, atoms[lastatom].type, angle);
        }
    }

    printf("Chain data:\n");
    for (i = 0; i < nchains; i++)
    {
        printf("Chain %d type:%d firstres:%d lastres:%d\n", i, chains[i].type, chains[i].firstres, chains[i].lastres);
    }

    return NO_ERROR;
}

