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
 * This file contains the definition of various functions that  are used to set molecule data when the chain type 
 * specifically corresponds to proteins. It builds atom connectivity (i.e.: bonds) assuming that the residues are one among
 * the 20 standard amino-acids.
 *
 */

#ifndef __PROTEIN_H__
#define __PROTEIN_H__

#include "datatypes.h"

#define char_is_digit(c)    ((47 < c) && (c < 58))

/**
 * Divides an atom name into its component parts and put them into the characters *at, *id, *n1 and *n2. 
 * For instance, for a string like "HA1", *at is set 'H', *id to 'A', *n1 to '1' and *n2 to ' '.
 * @param atname const char*
 * @param at char*
 * @param id char*
 * @param n1 char*
 * @param n2 char*
 * @return ErrorCode
 */
ErrorCode get_atom_name_parts(const char *atname, char *at, char *id, char *n1, char *n2);

/**
 * Converts the atom amino-acid into a type code value.
 * @param type TyepCode*
 * @param aaname const char*
 * @return ErrorCode
 */
ErrorCode get_aminoacid_code(TypeCode *type, const char *aaname);

/**
 * Converts the atom name into a type code value.
 * @param type TyepCode*
 * @param atname const char*
 * @return ErrorCode
 */
ErrorCode get_atom_code(TypeCode *type, const char *atname);
/**
 * Converts the type code value of an atom into a pdb name.
 * @param type TyepCode*
 * @param aaname char*
 * @return ErrorCode
 */
ErrorCode get_aminoacid_name(const TypeCode type,char *aaname);
/**
 * Converts the type code value into a atom name.
 * @param type TyepCode*
 * @param atname char*
 * @return ErrorCode
 */
ErrorCode get_atom_name(TypeCode type,char *atname);
/**
 * Returns in *atomidx the atom index of the atom in residue res with atomtype.
 * @param result IndexValue*
 * @param atomtype TypeCode
 * @param atoms struct AtomData*
 * @param residues struct ResidueData*
 * @param IndexValue res
 * @return ErrorCode
 */
ErrorCode get_atom_in_residue(IndexValue *atomidx, TypeCode atomtype, struct ResidueData *residues, IndexValue res, struct AtomData *atoms);

/**
 * Returns in *type the bond type corresponding to bond index bondidx.
 * @param type TypeCode*
 * @param bondidx IndexValue
 * @return ErrorCode
 */
ErrorCode get_bond_type(TypeCode *type, IndexValue bondidx);

/**
 * Sets the fields for bonds[bond].
 * @param bonds struct BondData*
 * @param bond IndexValue
 * @param atom0 IndexValue
 * @param atom1 IndexValue
 * @param atom2 IndexValue
 * @param atom3 IndexValue
 * @param lastatom IndexValue
 * @return ErrorCode
 */
ErrorCode set_valid_bond(struct BondData *bonds, IndexValue bond, IndexValue atom0, IndexValue atom1, IndexValue atom2, IndexValue atom3, IndexValue lastatom);

/**
 * Creates the main-chain bonds for residue res. It takes into account border conditions (residues located at the beginning or at the end of a chain).
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_mc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type ALA.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_ala_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type ALA. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_ala_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type ARG.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_arg_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type ARG. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_arg_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type ASN.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_asn_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type ASN. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_asn_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type ASP.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_asp_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type ASP. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_asp_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type CYS.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_cys_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type CYS. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_cys_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type GLN.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_gln_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type GLN. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_gln_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type GLU.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_glu_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type GLU. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_glu_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type GLY.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_gly_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type GLY. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_gly_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type HIS.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_his_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type HIS. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_his_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type ILE.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_ile_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type ILE. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_ile_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type LEU.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_leu_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type LEU. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_leu_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type LYS.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_lys_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type LYS. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_lys_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type MET.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_met_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type MET. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_met_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type PHE.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_phe_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type PHE. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_phe_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type PRO.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_pro_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type PRO. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_pro_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type SER.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_ser_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type SER. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_ser_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type THR.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_thr_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type THR. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_thr_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type TRP.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_trp_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type TRP. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_trp_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type TYR.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_tyr_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type TYR. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_tyr_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type VAL.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_val_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type VAL. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_val_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Creates the side-chain bonds for residue res, assuming that is of type UNK.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @return ErrorCode
 */
ErrorCode create_unk_sc_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains);
/**
 * Creates the bonds for residue res, assuming that is of type UNK. If createsc is FALSE, it only creates the main-chain bonds.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param atoms struct AtomData*
 * @param chains struct ChainData*
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_unk_bonds(struct ResidueData *residues, IndexValue res, struct AtomData *atoms, struct ChainData *chains, BoolValue createsc);

/**
 * Sets the fields of chains[chain] with the provided arguments.
 * @param chains struct ChainData*
 * @param chain IndexValue
 * @param firstres IndexValue
 * @param lastres IndexValue
 * @return ErrorCode
 */
ErrorCode set_protein_chain(struct ChainData *chains, IndexValue chain, IndexValue firstres, IndexValue lastres);

/**
 * Sets the fields of residues[res] with the provided arguments.
 * @param residues struct ResidueData*
 * @param res IndexValue
 * @param firstatom IndexValue
 * @param lastatom IndexValue
 * @param chain IndexValue
 * @param name const char*
 * @return ErrorCode
 */
ErrorCode set_protein_residue(struct ResidueData *residues, IndexValue res, IndexValue firstatom, IndexValue lastatom, IndexValue chain, const char *name);

/**
 * Sets the fields of atoms[atom] with the provided arguments.
 * @param residues struct AtomData*
 * @param atom IndexValue
 * @param x FloatValue
 * @param y FloatValue
 * @param z FloatValue
 * @param resIndexValue
 * @param chain IndexValue
 * @param name const char*
 * @return ErrorCode
 */
ErrorCode set_protein_atom(struct AtomData *atoms, IndexValue atom, FloatValue x, FloatValue y, FloatValue z, IndexValue res, IndexValue chain, const char *name);

/**
 * Creates the bonds for all the chains, assuming that they correspond to protein molecules. The createsc argument controls the creation of side-chain bonds.
 * @param atoms struct AtomData*
 * @param natoms IndexValue
 * @param residues struct ResidueData*
 * @param nres IndexValue
 * @param chains struct ChainData*
 * @param nchains IndexValue
 * @param createsc BoolValue
 * @return ErrorCode
 */
ErrorCode create_protein_bonds(struct AtomData *atoms, IndexValue natoms, struct ResidueData *residues, IndexValue nres, struct ChainData *chains, IndexValue nchains, BoolValue createsc);

#endif
