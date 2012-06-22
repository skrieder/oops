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

#ifndef __PROTCONSTANTS_H__
#define __PROTCONSTANTS_H__

#include "base.h"

/**
 * Amino-acid types.
 */
/*
static const TypeCode ALA = 0;
static const TypeCode ARG = 1;
static const TypeCode ASN = 2;
static const TypeCode ASP = 3;
static const TypeCode CYS = 4;
static const TypeCode GLN = 5;
static const TypeCode GLU = 6;
static const TypeCode GLY = 7;
static const TypeCode HIS = 8;
static const TypeCode ILE = 9;
static const TypeCode LEU = 10;
static const TypeCode LYS = 11;
static const TypeCode MET = 12;
static const TypeCode PHE = 13;
static const TypeCode PRO = 14;
static const TypeCode SER = 15;
static const TypeCode THR = 16;
static const TypeCode TRP = 17;
static const TypeCode TYR = 18;
static const TypeCode VAL = 19;
static const TypeCode UNK = 20;
*/

//comment out the above and uncomment this to use 1-numbered amino acids
static const TypeCode ALA =  1;
static const TypeCode ARG =  2;
static const TypeCode ASN =  3;
static const TypeCode ASP =  4;
static const TypeCode CYS =  5;
static const TypeCode GLN =  6;
static const TypeCode GLU =  7;
static const TypeCode GLY =  8;
static const TypeCode HIS =  9;
static const TypeCode ILE = 10;
static const TypeCode LEU = 11;
static const TypeCode LYS = 12;
static const TypeCode MET = 13;
static const TypeCode PHE = 14;
static const TypeCode PRO = 15;
static const TypeCode SER = 16;
static const TypeCode THR = 17;
static const TypeCode TRP = 18;
static const TypeCode TYR = 19;
static const TypeCode VAL = 20;
static const TypeCode UNK = 21;

/**
 * Atom types in a molecule chain.
 */
static const TypeCode MC_N   =  1;
static const TypeCode MC_CA  =  2;
static const TypeCode MC_C   =  3;
static const TypeCode MC_O   =  4;
static const TypeCode MC_HN  =  5;
static const TypeCode MC_HA1 =  6;
static const TypeCode MC_CAP =  7;

static const TypeCode SC_CB  =  8;
static const TypeCode SC_HA2 =  9;

static const TypeCode SC_CG  = 10;
static const TypeCode SC_CG1 = 11;
static const TypeCode SC_CG2 = 12;
static const TypeCode SC_OG  = 13;
static const TypeCode SC_OG1 = 14;
static const TypeCode SC_SG  = 15;

static const TypeCode SC_CD  = 16;
static const TypeCode SC_CD1 = 17;
static const TypeCode SC_CD2 = 18;
static const TypeCode SC_OD  = 19;
static const TypeCode SC_OD1 = 20;
static const TypeCode SC_OD2 = 21;
static const TypeCode SC_ND  = 22;
static const TypeCode SC_ND1 = 23;
static const TypeCode SC_ND2 = 24;
static const TypeCode SC_SD  = 25;

static const TypeCode SC_CE  = 26;
static const TypeCode SC_CE1 = 27;
static const TypeCode SC_CE2 = 28;
static const TypeCode SC_CE3 = 29;
static const TypeCode SC_OE  = 30;
static const TypeCode SC_OE1 = 31;
static const TypeCode SC_OE2 = 32;
static const TypeCode SC_NE  = 33;
static const TypeCode SC_NE1 = 34;
static const TypeCode SC_NE2 = 35;

static const TypeCode SC_CZ  = 36;
static const TypeCode SC_CZ2 = 37;
static const TypeCode SC_CZ3 = 38;
static const TypeCode SC_NZ  = 39;

static const TypeCode SC_CH2 = 40;
static const TypeCode SC_OH  = 41;
static const TypeCode SC_NH  = 42;
static const TypeCode SC_NH1 = 43; 
static const TypeCode SC_NH2 = 44;

static const TypeCode SC_H   = 45;
static const TypeCode SC_X   = 46;

/**
 * Bond types in a protein. For an abuse of notation, bonds are named after the torsional dihedral defined by them.
 */
static const TypeCode PHI_BOND       =  1;
static const TypeCode PSI_BOND       =  2;
static const TypeCode OMEGA_BOND     =  3;
static const TypeCode NH_BOND        =  4;
static const TypeCode CO_BOND        =  5;
static const TypeCode CAH_BOND       =  6;
static const TypeCode CHI1_BOND      =  7;
static const TypeCode CHI2_BOND      =  8;
static const TypeCode CHI3_BOND      =  9;
static const TypeCode CHI4_BOND      = 10;
static const TypeCode CHI5_BOND      = 11;
static const TypeCode UNDEFINED_BOND = 12;

/**
 * Indexes for protein bonds. For an abuse of notation, bonds are named after the torsional dihedral defined by them.
 */
static const IndexValue PHI   = 0;
static const IndexValue PSI   = 1;
static const IndexValue OMEGA = 2;
static const IndexValue CHI1  = 3;
static const IndexValue CHI2  = 4;
static const IndexValue CHI3  = 5;
static const IndexValue CHI4  = 6;
static const IndexValue CHI5  = 7;

#endif
