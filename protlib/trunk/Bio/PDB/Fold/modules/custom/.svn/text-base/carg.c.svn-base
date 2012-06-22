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
 * Sample custom module.
 *
 */

#include <stdio.h>
#include <math.h>
#include "../../core/memory.h"
#include "../../core/errorhandler.h"
#include "carg.h"

// Local copy of system variables.
struct AtomData *atoms;          // Array of atom data.
FloatValue *coords;              // Array of atom coordinates.
FloatValue *distances;           // Array of atom pair-wise distances.
FloatValue *invdistances;        // Array of atom pair-wise 1/distances.
BoolValue *atommask;             // Array of atom masking values.
BoolValue *distmask;             // Array of distance masking values.  
IndexValue natoms;               // Total number of atoms.

struct ResidueData *residues;    // Array of residue data.
IndexValue nres;                 // Total number of residues.
 
struct ChainData *chains;        // Array of chain data.
IndexValue nchains;              // Total number of chains.

IndexValue *cavec;

ErrorCode set_system_for_carg(struct AtomData *_atoms, FloatValue *_coords, 
                                                          FloatValue *_distances, FloatValue *_invdistances, 
                                                          BoolValue *_atommask, BoolValue *_distmask, IndexValue _natoms, 
                                 struct ResidueData *_residues, IndexValue _nres, 
                                 struct ChainData *_chains, IndexValue _nchains)
{
    atoms = _atoms;
    coords = _coords;
    distances = _distances;
    invdistances = _invdistances;
    atommask = _atommask;
    distmask = _distmask;
    natoms = _natoms;
    residues = _residues;
    nres = _nres;
    chains = _chains;
    nchains = _nchains;

    return NO_ERROR;
}

ErrorCode init_carg(const char* cfg)
{
    IndexValue i,count=0;
    printf("        Initializing Radius of gyration calculator...\n");
    create_index_array(&cavec,natoms);
    char str[16],str2[16];

    for(i=0;i<natoms;i++){
        if(atoms[i].type==MC_CA){
            cavec[count]=i;
            count++;
        }
    }
    sprintf(str," %i ",nres);
    sprintf(str2," %i ",count);
    if (count!=nres) print_error(MISMATCH,"carg_init",5,"Warning: number of residues (",str,") and number of c-alphas (",str2,") are not the same");

    printf("        Done.\n");
    return NO_ERROR;
}

ErrorCode finish_carg(void)
{
    delete_index_array(&cavec,natoms);
    return NO_ERROR;
}

ErrorCode pre_carg(void)
{
    return NO_ERROR;
}

ErrorCode post_carg(void)
{
    return NO_ERROR;
}

ErrorCode calc_cargb(BoolValue *carg, IndexValue n)
{
    return NO_ERROR;
}

ErrorCode calc_cargi(IntValue *carg, IndexValue n)
{
    return NO_ERROR;
}

ErrorCode calc_cargf(FloatValue *carg, IndexValue n)
{
    IndexValue i,j;
    FloatValue dist;

    FloatValue rad2=0;
    for(i=0;i<nres;i++)
        for(j=i+1;j<nres;j++){
            dist=distances[distance_index(cavec[i],cavec[j],natoms)];
            rad2+=(dist*dist); //note this is wasteful since in the distance matrix calculation it is already taking a square root
        }
    *carg = sqrt(rad2)/nres;

    return NO_ERROR;
}
