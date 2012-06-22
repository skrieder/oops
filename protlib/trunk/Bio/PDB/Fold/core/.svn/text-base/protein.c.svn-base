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
 * This file contains the implementation of various functions that  are used to set molecule data when the chain type 
 * specifically corresponds to proteins. It builds atom connectivity (i.e.: bonds) assuming that the residues are one among
 * the 20 standard amino-acids.
 *
 */

#include <string.h>
#include "memory.h"
#include "protein.h"
#include "errorhandler.h"

ErrorCode get_atom_name_parts(const char *atname, char *at, char *id, char *n1, char *n2)
{
    char c;
    IndexValue i;
    BoolValue prevletter;

    *at = *id = ' ';
    *n1 = *n2 = ' ';
    
    c = atname[0];
    i = 0;
    prevletter = FALSE;
    while (c != '\0')
    {
        if (char_is_digit(c) || (c == 'x'))
        {
            if (prevletter) *n1 = c; else *n2 = c;
            prevletter = FALSE;
        }
        else
        {
            if (prevletter) *id = c; else *at = c;
            prevletter = TRUE;
        }
        i++;
        c = atname[i];
    }

    if ((*n1 == ' ') && (*n2 != ' ')) { *n1 = *n2; *n2 = ' '; }

    return NO_ERROR;
}

ErrorCode get_aminoacid_code(TypeCode *type, const char *aaname)
{
    *type = UNK;
    if (strcmp(aaname, "ALA") == 0) *type = ALA;
    if (strcmp(aaname, "ARG") == 0) *type = ARG;
    if (strcmp(aaname, "ASN") == 0) *type = ASN;
    if (strcmp(aaname, "ASP") == 0) *type = ASP;
    if (strcmp(aaname, "CYS") == 0) *type = CYS;
    if (strcmp(aaname, "GLN") == 0) *type = GLN;
    if (strcmp(aaname, "GLU") == 0) *type = GLU;
    if (strcmp(aaname, "GLY") == 0) *type = GLY;
    if (strcmp(aaname, "HIS") == 0) *type = HIS;
    if (strcmp(aaname, "ILE") == 0) *type = ILE;
    if (strcmp(aaname, "LEU") == 0) *type = LEU;
    if (strcmp(aaname, "LYS") == 0) *type = LYS;
    if (strcmp(aaname, "MET") == 0) *type = MET;
    if (strcmp(aaname, "PHE") == 0) *type = PHE;
    if (strcmp(aaname, "PRO") == 0) *type = PRO;
    if (strcmp(aaname, "SER") == 0) *type = SER;
    if (strcmp(aaname, "THR") == 0) *type = THR;
    if (strcmp(aaname, "TRP") == 0) *type = TRP;
    if (strcmp(aaname, "TYR") == 0) *type = TYR;
    if (strcmp(aaname, "VAL") == 0) *type = VAL;
    return NO_ERROR;
}


ErrorCode get_aminoacid_name(const TypeCode type, char *aaname)
{
    if (type == ALA) strcpy(aaname, "ALA");
    if (type == ARG) strcpy(aaname, "ARG");
    if (type == ASN) strcpy(aaname, "ASN");
    if (type == ASP) strcpy(aaname, "ASP");
    if (type == CYS) strcpy(aaname, "CYS");
    if (type == GLN) strcpy(aaname, "GLN");
    if (type == GLU) strcpy(aaname, "GLU");
    if (type == GLY) strcpy(aaname, "GLY");
    if (type == HIS) strcpy(aaname, "HIS");
    if (type == ILE) strcpy(aaname, "ILE");
    if (type == LEU) strcpy(aaname, "LEU");
    if (type == LYS) strcpy(aaname, "LYS");
    if (type == MET) strcpy(aaname, "MET");
    if (type == PHE) strcpy(aaname, "PHE");
    if (type == PRO) strcpy(aaname, "PRO");
    if (type == SER) strcpy(aaname, "SER");
    if (type == THR) strcpy(aaname, "THR");
    if (type == TRP) strcpy(aaname, "TRP");
    if (type == TYR) strcpy(aaname, "TYR");
    if (type == VAL) strcpy(aaname, "VAL");
    return NO_ERROR;
}

ErrorCode get_atom_name(TypeCode type, char *atname)
{
    switch(type)
    {  // c numbers characters in octal!
        case '\1':  // MC_N
	    strcpy(atname,"N");
	    break;
        case '\2':  // MC_CA
	    strcpy(atname,"CA");
	    break;
	case '\3':  // MC_C
	    strcpy(atname,"C");
	    break;
        case '\4':  // MC_O
	    strcpy(atname,"O");
	    break;
        case '\5':  // MC_HN
	    strcpy(atname,  "HN");
	    break;
	case '\6':  // MC_HA1
	    strcpy(atname,  "HA");
	    break;
        case '\7':  // MC_CAP
	    strcpy(atname,  "OXT");
	    break;
        case '\10':  // SC_CB
	    strcpy(atname,  "CB");
	    break;
	case '\11':  // SC_HA2
	    strcpy(atname,  "HA2");
	    break;
        case '\12': // SC_CG
	    strcpy(atname,  "CG");
	    break;
        case '\13': // SC_CG1
	    strcpy(atname,  "CG1");
	    break;
	case '\14': // SC_CG2
	    strcpy(atname,  "CG2");
	    break;
        case '\15': // SC_OG
	    strcpy(atname,  "OG");
	    break;
        case '\16': // SC_OG1
	    strcpy(atname,  "OG1");
	    break;
	case '\17': // SC_SG
	    strcpy(atname,  "SG");
	    break;
        case '\20': // SC_CD
	    strcpy(atname,  "CD");
	    break;
        case '\21': // SC_CD1
	    strcpy(atname,  "CD1");
	    break;
	case '\22': // SC_CD2
	    strcpy(atname,  "CD2");
	    break;
	case '\23': // SC_OD
	    strcpy(atname,  "OD");
	    break;
        case '\24': // SC_OD1
	    strcpy(atname,  "OD1");
	    break;
        case '\25': // SC_OD2
	    strcpy(atname,  "OD2");
	    break;
	case '\26': // SC_ND
	    strcpy(atname,  "ND");
	    break;
        case '\27': // SC_ND1
	    strcpy(atname,  "ND1");
	    break;
        case '\30': // SC_ND2
	    strcpy(atname,  "ND2");
	    break;
	case '\31': // SC_SD
	    strcpy(atname,  "SD");
	    break;
        case '\32': // SC_CE
	    strcpy(atname,  "CE");
	    break;
        case '\33': // SC_CE1
	    strcpy(atname,  "CE1");
	    break;
	case '\34': // SC_CE2
	    strcpy(atname,  "CE2");
	    break;
        case '\35': // SC_CE3
	    strcpy(atname,  "CE3");
	    break;
        case '\36': // SC_OE
	    strcpy(atname,  "OE");
	    break;
	case '\37': // SC_OE1
	    strcpy(atname,  "OE1");
	    break;
        case '\40': // SC_OE2
	    strcpy(atname,  "OE2");
	    break;
        case '\41': // SC_NE
	    strcpy(atname,  "NE");
	    break;
	case '\42': // SC_NE1
	    strcpy(atname,  "NE1");
	    break;
        case '\43': // SC_NE2
	    strcpy(atname,  "NE2");
	    break;
        case '\44': // SC_CZ
	    strcpy(atname,  "CZ");
	    break;
	case '\45': // SC_CZ2
	    strcpy(atname,  "CZ2");
	    break;
        case '\46': // SC_CZ3
	    strcpy(atname,  "CZ3");
	    break;
        case '\47': // SC_NZ
	    strcpy(atname,  "NZ");
	    break;
	case '\50': // SC_CH2
	    strcpy(atname,  "CH2");
	    break;
        case '\51': // SC_OH
	    strcpy(atname,  "OH");
	    break;
        case '\52': // SC_NH
	    strcpy(atname,  "NH");
	    break;
	case '\53': // SC_NH1
	    strcpy(atname,  "NH1");
	    break;
	case '\54': // SC_NH2
	    strcpy(atname,  "NH2");
	    break;
	case '\55': // SC_H
	    strcpy(atname,  "H");
	    break;
	case '\56': // SC_X
	    strcpy(atname,  "X");
	    break;
    }
    return NO_ERROR;
}

ErrorCode get_atom_code(TypeCode *type, const char *atname)
{
    *type = SC_X;

    // Main-chain atoms.
    if (strcmp(atname, "N"  ) == 0) *type = MC_N;
    if (strcmp(atname, "CA" ) == 0) *type = MC_CA;
    if (strcmp(atname, "C"  ) == 0) *type = MC_C;
    if (strcmp(atname, "O"  ) == 0) *type = MC_O;
    if (strcmp(atname, "H"  ) == 0 || 
        strcmp(atname, "HN" ) == 0) *type = MC_HN;
    if (strcmp(atname, "HA" ) == 0 || 
        strcmp(atname, "HA1") == 0 ||
        strcmp(atname, "1HA") == 0) *type = MC_HA1;

    // Capping atoms.
    if (strcmp(atname, "HT1") == 0 || 
        strcmp(atname, "HT2") == 0 ||
        strcmp(atname, "HT3") == 0 || 
        strcmp(atname, "1H" ) == 0 ||
        strcmp(atname, "2H" ) == 0 || 
        strcmp(atname, "3H" ) == 0 ||
        strcmp(atname, "OXT") == 0) *type = MC_CAP;

    // Main-chain atoms.
    if (strcmp(atname, "CB" ) == 0) *type = SC_CB;
    if (strcmp(atname, "HA2") == 0 || 
        strcmp(atname, "2HA") == 0) *type = SC_HA2;

    // Side-chain atoms.
    if (strcmp(atname, "CG" ) == 0) *type = SC_CG;
    if (strcmp(atname, "CG1") == 0) *type = SC_CG1;
    if (strcmp(atname, "CG2") == 0) *type = SC_CG2;
    if (strcmp(atname, "OG" ) == 0) *type = SC_OG;
    if (strcmp(atname, "OG1") == 0) *type = SC_OG1;
    if (strcmp(atname, "SG" ) == 0) *type = SC_SG;

    if (strcmp(atname, "CD" ) == 0) *type = SC_CD;
    if (strcmp(atname, "CD1") == 0) *type = SC_CD1;
    if (strcmp(atname, "CD2") == 0) *type = SC_CD2;
    if (strcmp(atname, "OD" ) == 0) *type = SC_OD;
    if (strcmp(atname, "OD1") == 0) *type = SC_OD1;
    if (strcmp(atname, "OD2") == 0) *type = SC_OD2;
    if (strcmp(atname, "ND" ) == 0) *type = SC_ND;
    if (strcmp(atname, "ND1") == 0) *type = SC_ND1;
    if (strcmp(atname, "ND2") == 0) *type = SC_ND2;
    if (strcmp(atname, "SD" ) == 0) *type = SC_SD;

    if (strcmp(atname, "CE" ) == 0) *type = SC_CE;
    if (strcmp(atname, "CE1") == 0) *type = SC_CE1;
    if (strcmp(atname, "CE2") == 0) *type = SC_CE2;
    if (strcmp(atname, "CE3") == 0) *type = SC_CE3;
    if (strcmp(atname, "OE" ) == 0) *type = SC_OE;
    if (strcmp(atname, "OE1") == 0) *type = SC_OE1;
    if (strcmp(atname, "OE2") == 0) *type = SC_OE2;
    if (strcmp(atname, "NE" ) == 0) *type = SC_NE;
    if (strcmp(atname, "NE1") == 0) *type = SC_NE1;
    if (strcmp(atname, "NE2") == 0) *type = SC_NE2;

    if (strcmp(atname, "CZ" ) == 0) *type = SC_CZ;
    if (strcmp(atname, "CZ2") == 0) *type = SC_CZ2;
    if (strcmp(atname, "CZ3") == 0) *type = SC_CZ3;
    if (strcmp(atname, "NZ" ) == 0) *type = SC_NZ;

    if (strcmp(atname, "CH2") == 0) *type = SC_CH2;
    if (strcmp(atname, "OH" ) == 0) *type = SC_OH;
    if (strcmp(atname, "NH" ) == 0) *type = SC_NH;
    if (strcmp(atname, "NH1") == 0) *type = SC_NH1;
    if (strcmp(atname, "NH2") == 0) *type = SC_NH2;
    
    if (*type == SC_X)
    {
        char at, id, n1, n2;
        get_atom_name_parts(atname, &at, &id, &n1, &n2);
        if (at == 'H') *type = SC_H;
    }

    return NO_ERROR;
}

ErrorCode get_atom_in_residue(IndexValue *atomidx, TypeCode atomtype, struct ResidueData *residues, IndexValue res, struct AtomData *atoms)
{
    IndexValue firstatom = residues[res].firstatom;
    IndexValue lastatom = residues[res].lastatom;
    
    IndexValue atom;
    for (atom = firstatom; atom <= lastatom; atom++)
        if (atomtype == atoms[atom].type)
        {
            *atomidx = atom;
            return NO_ERROR;
        }
    
    return ITEM_NOT_FOUND_ERROR;
}

ErrorCode get_bond_type(TypeCode *type, IndexValue bondidx)
{
    *type = UNDEFINED_BOND;
    if (bondidx == PHI  ) *type = PHI_BOND;
    if (bondidx == PSI  ) *type = PSI_BOND;
    if (bondidx == OMEGA) *type = OMEGA_BOND;
    if (bondidx == CHI1 ) *type = CHI1_BOND;
    if (bondidx == CHI2 ) *type = CHI2_BOND;
    if (bondidx == CHI3 ) *type = CHI3_BOND;
    if (bondidx == CHI4 ) *type = CHI4_BOND;
    if (bondidx == CHI5 ) *type = CHI5_BOND;

    return NO_ERROR;
}

ErrorCode set_valid_bond(struct BondData *bonds, IndexValue bond, IndexValue atom0, IndexValue atom1, IndexValue atom2, IndexValue atom3, IndexValue lastatom)
{
    get_bond_type(&(bonds[bond].type), bond);

    bonds[bond].atom0 = atom0;
    bonds[bond].atom1 = atom1;
    bonds[bond].atom2 = atom2;
    bonds[bond].atom3 = atom3;
 
    bonds[bond].lastatom = lastatom;

    return NO_ERROR;
}

ErrorCode create_mc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    IndexValue firstchainres = chains[residues[res].chain].firstres;
    IndexValue lastchainres = chains[residues[res].chain].lastres;
    IndexValue lastchainatom = residues[lastchainres].lastatom;

    BoolValue firstres = TRUE;
    IndexValue prevres = 0;
    if (firstchainres < res)
    {
        prevres = res - 1;
        firstres = FALSE;
    }

    BoolValue lastres = TRUE;
    IndexValue nextres = 0;
    if (res < lastchainres)
    {
        nextres = res + 1;
        lastres = FALSE;
    }

    if (!firstres)
    {
        error0 = get_atom_in_residue(&atom0, MC_C, residues, prevres, atoms);
        error1 = get_atom_in_residue(&atom1, MC_N, residues, res, atoms);
        error2 = get_atom_in_residue(&atom2, MC_CA, residues, res, atoms);
        error3 = get_atom_in_residue(&atom3, MC_C, residues, res, atoms);
        if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
        {
            set_valid_bond(residues[res].bonds, PHI, atom0, atom1, atom2, atom3, lastchainatom);
        }
    }

    if (!lastres)
    {
        error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
        error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
        error2 = get_atom_in_residue(&atom2, MC_C, residues, res, atoms);
        error3 = get_atom_in_residue(&atom3, MC_N, residues, nextres, atoms);
        if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
        {
            set_valid_bond(residues[res].bonds, PSI, atom0, atom1, atom2, atom3, lastchainatom);
        }
    }

    if (!lastres)
    {
        error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
        error1 = get_atom_in_residue(&atom1, MC_C, residues, res, atoms);
        error2 = get_atom_in_residue(&atom2, MC_N, residues, nextres, atoms);
        error3 = get_atom_in_residue(&atom3, MC_CA, residues, nextres, atoms);
        if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
        {
            set_valid_bond(residues[res].bonds, OMEGA, atom0, atom1, atom2, atom3, lastchainatom);
        }
    }

    return NO_ERROR;
}

ErrorCode create_ala_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    return NO_ERROR;
}

ErrorCode create_ala_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 0;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_ala_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_arg_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_NH2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, SC_CB, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CG, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CD, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_NE, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI3, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, SC_CG, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CD, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_NE, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CZ, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI4, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, SC_CD, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_NE, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_NZ, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_NH1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI5, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_arg_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 5;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_arg_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_asn_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_ND2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_OD1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_asn_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 2;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_asn_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_asp_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_OD2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_OD1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_asp_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 2;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_asp_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_cys_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_SG, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_SG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_cys_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 1;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_cys_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_gln_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_NE2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, SC_CB, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CG, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CD, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_OE1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI3, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_gln_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 3;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_gln_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_glu_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_OE2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, SC_CB, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CG, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CD, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_OE1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI3, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_glu_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 3;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_glu_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_gly_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    return NO_ERROR;
}

ErrorCode create_gly_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 0;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_gly_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_his_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_NE2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_ND1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_his_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 2;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_his_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_ile_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_CD1, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG1, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_ile_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 2;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_ile_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_leu_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_CD2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_leu_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 2;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_leu_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_lys_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_NZ, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, SC_CB, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CG, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CD, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CE, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI3, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, SC_CG, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CD, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CE, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_NZ, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI4, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_lys_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 4;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_lys_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_met_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_CE, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_SD, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, SC_CB, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CG, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_SD, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CE, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI3, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_met_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 3;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_met_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_phe_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_CZ, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_phe_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 2;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_phe_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_pro_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    return NO_ERROR;
}

ErrorCode create_pro_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 0;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_pro_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_ser_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_OG, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_OG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_ser_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 1;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_ser_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_thr_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_CG2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_OG1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_thr_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 1;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_thr_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}


ErrorCode create_trp_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_CH2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_trp_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 2;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_trp_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_tyr_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_OH, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    error0 = get_atom_in_residue(&atom0, MC_CA, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, SC_CB, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CG, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CD1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI2, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_tyr_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 2;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_tyr_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_val_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    ErrorCode error0, error1, error2, error3;
    IndexValue atom0, atom1, atom2, atom3; 

    ErrorCode error;

    IndexValue lastatom;
    error = get_atom_in_residue(&lastatom, SC_CG2, residues, res, atoms);
    if (error != NO_ERROR) return error;

    error0 = get_atom_in_residue(&atom0, MC_N, residues, res, atoms);
    error1 = get_atom_in_residue(&atom1, MC_CA, residues, res, atoms);
    error2 = get_atom_in_residue(&atom2, SC_CB, residues, res, atoms);
    error3 = get_atom_in_residue(&atom3, SC_CG1, residues, res, atoms);
    if (error0 == NO_ERROR && error1 == NO_ERROR && error2 == NO_ERROR && error3 == NO_ERROR)
    {
        set_valid_bond(residues[res].bonds, CHI1, atom0, atom1, atom2, atom3, lastatom);
    }

    return NO_ERROR;
}

ErrorCode create_val_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 1;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_val_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode create_unk_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains)
{
    return NO_ERROR;
}

ErrorCode create_unk_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc)
{
    ErrorCode error;

    residues[res].nbonds = 3;
    if (createsc) residues[res].nbonds += 0;
    error = create_bonds_array(&(residues[res].bonds), residues[res].nbonds);
    if (error != NO_ERROR) return error;

    error = create_mc_bonds(residues, res, atoms, chains);
    if (error != NO_ERROR) return error;

    if (createsc)
    {
        error = create_unk_sc_bonds(residues, res, atoms, chains);
        if (error != NO_ERROR) return error;
    }

    return NO_ERROR;
}

ErrorCode set_protein_chain(struct ChainData *chains, IndexValue chain, IndexValue firstres, IndexValue lastres)
{
    chains[chain].type = PROTEIN;

    chains[chain].firstres = firstres;
    chains[chain].lastres = lastres;

    return NO_ERROR;
}

ErrorCode set_protein_residue(struct ResidueData *residues, IndexValue res, IndexValue firstatom, IndexValue lastatom, IndexValue chain, const char *name)
{
    get_aminoacid_code(&(residues[res].type), name);

    residues[res].firstatom = firstatom;
    residues[res].lastatom = lastatom;

    residues[res].chain = chain;

    return NO_ERROR;
}
 
ErrorCode set_protein_atom(struct AtomData *atoms, IndexValue atom, FloatValue x, FloatValue y, FloatValue z, IndexValue res, IndexValue chain, const char *name)
{
    get_atom_code(&(atoms[atom].type), name);

    *(atoms[atom].x) = x;
    *(atoms[atom].y) = y;
    *(atoms[atom].z) = z;

    atoms[atom].res = res;
    atoms[atom].chain = chain;

    return NO_ERROR;
}

ErrorCode create_protein_bonds(struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, IndexValue nres, struct ChainData *chains, IndexValue nchains, BoolValue createsc)
{
    ErrorCode error0, error;
    IndexValue res;
   
    error0 = NO_ERROR;
    for (res = 0; res < nres; res++)
    {
        if (residues[res].type == ALA) 
        {
            error = create_ala_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_ala_bonds", 0);
        }
        else if (residues[res].type == ARG) 
        {
            error = create_arg_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_arg_bonds", 0);
        }
        else if (residues[res].type == ASN) 
        {
            error = create_asn_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_asn_bonds", 0);
        }
        else if (residues[res].type == ASP) 
        {
            error = create_asp_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_asp_bonds", 0);
        }
        else if (residues[res].type == CYS) 
        {
            error = create_cys_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_cys_bonds", 0);
        }
        else if (residues[res].type == GLN) 
        {
            error = create_gln_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_gln_bonds", 0);
        }
        else if (residues[res].type == GLU) 
        {
            error = create_glu_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_glu_bonds", 0);
        }
        else if (residues[res].type == GLY) 
        {
            error = create_gly_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_gly_bonds", 0);
        }
        else if (residues[res].type == HIS) 
        {
            error = create_his_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_his_bonds", 0);
        }
        else if (residues[res].type == ILE) 
        {
            error = create_ile_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_ile_bonds", 0);
        }
        else if (residues[res].type == LEU) 
        {
            error = create_leu_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_leu_bonds", 0);
        }
        else if (residues[res].type == LYS) 
        {
            error = create_lys_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_lys_bonds", 0);
        }
        else if (residues[res].type == MET) 
        {
            error = create_met_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_met_bonds", 0);
        }
        else if (residues[res].type == PHE) 
        {
            error = create_phe_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_phe_bonds", 0);
        }
        else if (residues[res].type == PRO) 
        {
            error = create_pro_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_pro_bonds", 0);
        }
        else if (residues[res].type == SER) 
        {
            error = create_ser_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_ser_bonds", 0);
        }
        else if (residues[res].type == THR) 
        {
            error = create_thr_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_thr_bonds", 0);
        }
        else if (residues[res].type == TRP) 
        {
            error = create_trp_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_trp_bonds", 0);
        }
        else if (residues[res].type == TYR) 
        {
            error = create_tyr_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_tyr_bonds", 0);
        }
        else if (residues[res].type == VAL) 
        {
            error = create_val_bonds(residues, res, atoms, chains, createsc);
            print_error(error, "create_val_bonds", 0);
        }
        else                                
        {
            error = create_unk_bonds(residues, res, atoms, chains, createsc); 
            print_error(error, "create_unk_bonds", 0);
        }
        
        if (error) error0 = BUILDING_BOND_ERROR;
    }
    return error0;
}
