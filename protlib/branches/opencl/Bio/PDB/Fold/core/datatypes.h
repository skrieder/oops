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
 * This file contains the declaration of basic data structure types to hold molecule information.
 * There is no assumption made here as to the specific type of molecules that will be represented
 * with these structures. There are only two assumptions made with regard of the structure fo the
 * data that will be stored in these structures: 
 * 1) the hierarchy of chain->residue->atom, although chain could be protein, rna, dna, etc. For protein 
 * chains, residues should be amino-acids, nucleic bases for rna/dna, etc.
 * 2) the ordering of the atoms of a branched molecule in the array of AtomData structures. The numbering
 * of the atoms within a given branch should be consecutive, and the atoms of the branches starting at a 
 * given atom should occur after this atom. In order words, if i is the array index of a branching atom,
 * and branches B0, B1, B2, etc, start at b, then the branches should be all the atoms within an index interval
 * [j00, j01] for B0, [j10, j11] for B1 and so on, Furthermore:
 * i, j00 = i + 1, ..., j01, j10 = j01 + 1, ..., j11...
 * In this way, torsional rotations can be defined and implemented in a very convenient way.
 *
 */

#ifndef __DATATYPES_H__
#define __DATATYPES_H__

#include "base.h"
#include "errorcodes.h"
#include "constants.h"

/**
 * Record of atom data.
 */
struct AtomData
{
    /**
            * Atom type (Carbon, Hydrogen, Oxygen, etc).
            */
    TypeCode type;

    /**
            * Pointers to the xyzw coordinates of the atom.
            */
    FloatValue *x;
    FloatValue *y;
    FloatValue *z;
    FloatValue *w;

    /**
            * Residue and chain where this atom resides.
            */
    IndexValue res;
    IndexValue chain;

    /**
            * Arrays to hold custom integer or float properties for each atom.
            * This arrays are always NULL when the atom record is initialized. 
            */
    IndexValue niprop;
    IntValue *iproperties;
    IndexValue nfprop;
    FloatValue *fproperties;
};

/**
 * Record of bond data. A bond is defined by four consecutive atoms in a chain, and 
 * have a torsional or dihedral angle associated to it. By rotating this angle, the atoms
 * that occur down the chain starting from index atom3 will change their position 
 * in space. 
 */
struct BondData
{
    /**
            * Bond type.
            */
    TypeCode type;

    /**
            * Indexes of the atoms defining the bond. It should be atom0 < atom1 < atom2 < atom3,
            * but not necessarily atom1 = atom0 + 1, etc.
            */
    IndexValue atom0; 
    IndexValue atom1; 
    IndexValue atom2;
    IndexValue atom3;

    /**
            * Index of the last atom in the chain that is affected by a torsional rotation at this bond.
            * The existence of such an index requires the atom ordering described at beginning of
            * this file, this is, i, j00 = i + 1, ..., j01, j10 = j01 + 1, ..., j11... for branches defined at i.
            */
    IndexValue lastatom;

    /**
            * Current value for the torsional angle of this bond.
            */
    FloatValue angle;
    /**
            * Used to store a new value for the torsional angle.
            */
    FloatValue newangle;
};

/**
 * Record of residue data.
 */
struct ResidueData
{
    /**
            * Residue type. This should represent amino-acids in the case of protein chains, nucleic acids
            * for rna, etc.
            */
    TypeCode type;

    /**
            * First and last atom of this residue.
            */
    IndexValue firstatom;
    IndexValue lastatom;

    /**
            * Chain containing this residue.
            */     
    IndexValue chain; 

    /**
            * Numer of bonds that start at atoms within this residue.
            */
    IndexValue nbonds;
    /**
            * List of bonds that start at atoms within this residue.
            */
    struct BondData *bonds;

    /**
            * Arrays to hold custom integer or float properties for each residue.
            * This arrays are always NULL when the residue record is initialized. 
            */
    IndexValue niprop;
    IntValue *iproperties;
    IndexValue nfprop;
    FloatValue *fproperties;
};

/**
 * Record of chain data.
 */
struct ChainData
{
    /**
            * Chain type.
            */
    TypeCode type;

    /**
            * First and last residue of this chain.
            */
    IndexValue firstres;
    IndexValue lastres;
};

// Some useful defines to access the components of an atom in a 1D array where the coordinates are stored one
// after another.
#define xcoord_index(atom)	(atom * 4 + 0)
#define ycoord_index(atom)	(atom * 4 + 1)
#define zcoord_index(atom)	(atom * 4 + 2)
#define wcoord_index(atom)	(atom * 4 + 3)

// Similarly, a convenience define to access the elements of a 2D  distance matrix stored as a 1D array.
#define distance_index(atom0, atom1, natoms)	(atom0 * natoms + atom1)

#endif