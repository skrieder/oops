# The C Protein Folding Library.
# Copyright (C) 2009 Andres Colubri.
# Contact: andres.colubri 'AT' gmail.com
#
# This library was written at the Institute of Biophysical Dynamics at the University of Chicago.
# Gordon Center for Integrated Science, W101. 929 East 57th Street, Chicago, IL 60637, USA.
# Homepage: http://ibd.uchicago.edu/
# 
# Permission to use, copy, modify, and distribute this software and its
# documentation with or without modifications and for any purpose and
# without fee is hereby granted, provided that any copyright notices
# appear in all copies and that both those copyright notices and this
# permission notice appear in supporting documentation, and that the
# names of the contributors or copyright holders not be used in
# advertising or publicity pertaining to distribution of the software
# without specific prior permission.
# 
# THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
# WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
# CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
# OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
# OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
# OR PERFORMANCE OF THIS SOFTWARE. 

# This file implements the interface objects and functions between Biopython and the C Protein Folding Library.
# Biopython should be in charge of handling the input and output of structures. At this point, only protein structures
# are supported.

import numpy
from types import StringType

from Bio.PDB import *
from fold import *
import build_from_seq
from distance_matrix import make_dij
import rmsd

PROTEIN = 1

# Atom types currently recognized by the FoldingSimulator object.
supported_atoms = {"N"  : 0, 
                   "HN" : 1, "H"  : 2, 
                   "CA" : 3, 
                   "CB" : 4,
                   "CG" : 5, "CG1": 6, "CG2": 7, "OG" : 8, "OG1": 9, "SG" :10,
                   "CD" :11, "CD1":12, "CD2":13, "OD" :14, "OD1":15, "OD2":16, "ND" :17, "ND1":18, "ND2":19, "SD" :20,
                   "CE" :21, "CE1":22, "CE2":23, "CE3":24, "OE" :25, "OE1":26, "OE2":27, "NE" :28, "NE1":29, "NE2":30,
                   "CZ" :31, "CZ2":32, "CZ3":33, "NZ" :34,
                   "CH2":35, "OH" :36, "NH" :37, "NH1":38, "NH2":39,
                   "HA1":40, "1HA":41, "HA" :42,
                   "HA2":43, "2HA":44,
                   "C"  :45, "O"  :46}

# Dictionary that specify the  order in which the atoms should be passed to the underlying C library in the case of not using full side-chains.
aa_nosc_atoms = {"N":0, "HN":1, "CA":2, "CB":3, "HA2":4, "HA1":5, "C":6, "O":7}

# These dictionaries specify the order in which the atoms should be passed to the underlying C library for each amino-acid. This ordering is important to properly
# construct the protein connectivity.
ala_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "HA1":4, "C"  :5, "O"  :6}
arg_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "CD" :5, "NE" :6, "CZ" :7, "NH1":8, "NH2":9, "HA1":10, "C"  :11,  "O" :12}
asn_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "OD1":5, "ND2":6, "HA1":7, "C"  :8, "O"  :9}
asp_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "OD1":5, "OD2":6, "HA1":7, "C"  :8, "O"  :9}
cys_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "SG" :4, "HA1":5, "C"  :6, "O"  :7}
gln_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "CD" :5, "OE1":6, "NE2":7, "HA1":8, "C"  :9, "O"  :10}
glu_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "CD" :5, "OE1":6, "OE2":7, "HA1":8, "C"  :9, "O"  :10}
gly_atoms = {"N":0, "HN":1, "CA":2, "HA2":3, "HA1":4, "C"  :5, "O"  :6}
his_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "ND1":5, "CD2":6, "CE1":7, "NE2":8, "HA1":9, "C"  :10, "O"  :11}
ile_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG1":4, "CD1":5, "CG2":6, "HA1":7, "C"  :8, "O"  :9}
leu_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "CD1":5, "CD2":6, "HA1":7, "C"  :8, "O"  :9}
lys_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "CD" :5, "CE" :6, "NZ" :7, "HA1":8, "C"  :9, "O"  :10}
met_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "SD" :5, "CE" :6, "HA1":7, "C"  :8, "O"  :9}
phe_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "CD1":5, "CD2":6, "CE1":7, "CE2":8, "CZ" :9, "HA1":10, "C"  :11,  "O" :12}
pro_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "HA1":5, "C"  :6, "O"  :7}
ser_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "OG" :4, "HA1":5, "C"  :6, "O"  :7}
thr_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "OG1":4, "CG2":5, "HA1":6, "C"  :7, "O"  :8}
trp_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "CD1":5, "CD2":6, "NE1":7, "CE2":8, "CE3":9, "CZ2":10, "CZ3":11, "CH2":12, "HA1":13, "C" :14,  "O" :15}
tyr_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG" :4, "CD1":5, "CD2":6, "CE1":7, "CE2":8, "CZ" :9, "OH" :10, "HA1":11, "C"  :12,  "O" :13}
val_atoms = {"N":0, "HN":1, "CA":2, "CB" :3, "CG1":4, "CG2":5, "HA1":6, "C"  :7, "O"  :8}

aa_atoms = [ala_atoms, arg_atoms, asn_atoms, asp_atoms, cys_atoms, gln_atoms, glu_atoms, gly_atoms, his_atoms, ile_atoms,
            leu_atoms, lys_atoms, met_atoms, phe_atoms, pro_atoms, ser_atoms, thr_atoms, trp_atoms, tyr_atoms, val_atoms]

# Aminoacid-to-index lookup table
aa_names = {"ALA": 0, "ARG": 1, "ASN": 2, "ASP": 3, "CYS": 4, "GLN": 5, "GLU": 6, "GLY": 7, "HIS": 8, "ILE": 9,
            "LEU":10, "LYS":11, "MET":12, "PHE":13, "PRO":14, "SER":15, "THR":16, "TRP":17, "TYR":18, "VAL":19}

#Default atom codes for protlib2
aa_def_types = {'MC_N': 1,'MC_CA': 2,'MC_C': 3,'MC_O': 4,'MC_HN': 5,'MC_HA1': 6,'MC_CAP': 7,'SC_CB': 8,'SC_HA2': 9,'SC_CG': 10,'SC_CG1': 11,'SC_CG2': 12,'SC_OG': 13,'SC_OG1': 14,'SC_SG': 15,'SC_CD': 16,'SC_CD1': 17,'SC_CD2': 18,'SC_OD': 19,'SC_OD1': 20,'SC_OD2': 21,'SC_ND': 22,'SC_ND1': 23,'SC_ND2': 24,'SC_SD': 25,'SC_CE': 26,'SC_CE1': 27,'SC_CE2': 28,'SC_CE3': 29,'SC_OE': 30,'SC_OE1': 31,'SC_OE2': 32,'SC_NE': 33,'SC_NE1': 34,'SC_NE2': 35,'SC_CZ': 36,'SC_CZ2': 37,'SC_CZ3': 38,'SC_NZ': 39,'SC_CH2': 40,'SC_OH': 41,'SC_NH': 42,'SC_NH1': 43,'SC_NH2': 44,'SC_H': 45,'SC_X': 46}

def seq_to_PDB(seqfile,format="fasta",datafile=None,returndict=False):
    """
    This function uses the submodule build_from_seq and techniques from LINUS to construct a three dimensional structure out of a protein sequence file (or open file handle).  Data file must be in the format of ribosome.dat from LINUS is automatically acquired if not specified. Return value is a dictionary of structures if multiple proteins were in the input and otherwise is a single PDBparser object. 

    The optional argument returndict can be set as True to return a dictionary regardless of the number of sequences.
    """
    protdict=build_from_seq.seq_to_PDB(seqfile,format=format,paramfile=datafile)
    if len(protdict.keys())>1 or returndict:
        return protdict
    else:
        return protdict[protdict.keys()[0]]

def aaname_to_index(s):
    """
    Three letter code to index.
    For example: ALA to 0.
    """
    return aa_names[s]

def normalize_hydrogen_name(name):
    """
    Convert the name of an hydrogen atom into a standard format.
    For example: H to HN.
    """
    if name == "H":
        return "HN"

    if name == "1HA" or name == "HA":
        return "HA1"

    if name == "2HA":
        return "HA2"

    return name

def is_supported_at(atom):
    """
    Return 1 if atom object/string is a supported atom.

    @param atom: a {Atom} object OR a atom name
    @type atom: {Atom} or string

    @type standard: boolean
    """
    if not type(atom) == StringType:
        atom = atom.get_name()
    atom = atom.upper()

    return supported_atoms.has_key(atom)

class FoldingSimulator:
    """
    An object that encapsulates the C Protein Folding Library and allows to run
    folding simulations from Biopython. A Bio.PDB.Structure object must be passed
    to the object with set_structure, and the atomic coordinates of this structure 
    are copied to the underlying C library. This object provides methods to set the 
    simulation driver. After the simulation is completed, the resulting coordinates 
    are copied back to the structure object.
    """
    def __init__(self):
        """
        Default constructor.
        """
        self.proteins = None

    def set_structure(self, structure):
        """
        Sets the Structure object to use as initial state for the simulation.
        """
        self.proteins = structure
        pl_set_driver_by_index(0)

    def get_drivers(self):
        """
        Returns a list with the name of the currently available driver modules.
        """
        ndrv = pl_get_num_drivers()
        list = []
        for n in range(0, ndrv):
            list.append(pl_get_driver_name(n))
        return list

    def get_driver_doc(self, drv):
        """
        Returns the doc string for the specified driver.
        """
        return pl_get_driver_doc(drv)

    def set_driver(self, drv):
        """
        Sets the driver drv as current. drv can be a string (the name of the driver) 
        or an integer number (the index of the driver).
        """
        if type(drv) == StringType:
            pl_set_driver_by_name(drv)
        else:
            pl_set_driver_by_index(drv)
    
    def run_simulation(self, cfg, use_model=0, use_sc=0):
        """
        Runs a folding simulation, using cfg as the configuration file for the current
        driver. The arguments use_model and use_sc allows to select which model in the 
        protein structure to use, and to enable/disable the use of full side-chains, 
        respectively. Non full-side-chains mode load the Beta Carbon.
        The resulting conformation after the end of the simulation is copied back to 
        the structure object.
        """
        model = self.proteins[use_model]

        # This dictionary maps the atom indexes in the underlying C library to the atom objects in the Biopython structure object.
        atom_dict = {}

        # These lists store the first and last atom for each residue.
        firstatom = []
        lastatom = []

        # These lists store the first and last residue for each chain.
        firstres = []
        lastres = []

        nchains = 0
        nres = 0
        natoms = 0
        chainidx = 0
        residx = 0
        atomidx = 0

        # Getting total number of chains, residues and atoms.
        print "CALCULATING TOTAL NUMBER OF ATOMS, RESIDUES AND CHAINS..."
        for chain in model:
            nchains = nchains + 1
            for res in chain:
                aaname = res.get_resname()
                if is_aa(aaname):
                    aa = aaname_to_index(aaname)
                    nres = nres + 1
                    for atom in res:
                        if is_supported_at(atom):
                            atname = normalize_hydrogen_name(atom.get_name())
                            if use_sc:
                                atdict = aa_atoms[aa]
                            else:
                                atdict = aa_nosc_atoms                            
                            if atdict.has_key(atname):
                                natoms = natoms + 1

        print "    ATOMS   : ", natoms
        print "    RESIDUES: ", nres
        print "    CHAINS  : ", nchains

        # Creating the data arrays in C.
        print "GETTING MEMORY FOR SIMULATION...",
        pl_allocate_memory(natoms, nres, nchains)
        print "DONE."

        # Copying data from Biopython structure object to the underlying C variables (atoms, residues and chains).
        print "COPYING DATA FROM BIOPYTHON...",
        for chain in model:
            firstres = firstres + [residx]
            for res in chain:
                aaname = res.get_resname()
                if is_aa(aaname):
                    aa = aaname_to_index(aaname)
                    if use_sc:
                        atdict = aa_atoms[aa]
                    else:
                        atdict = aa_nosc_atoms
                    nat = len(atdict)

                    # Getting all the supported atoms that are present in this residue.
                    present_atoms = [None] * nat
                    for atom in res:
                        if is_supported_at(atom):
                            atname = normalize_hydrogen_name(atom.get_name())
                            if atdict.has_key(atname):
                                at = atdict[atname]
                                present_atoms[at] = atom

                    # Copying atom data to C variables.
                    atomidx0 = atomidx
                    for at in range(0, nat):
                        if present_atoms[at] != None:
                            atom = present_atoms[at]
                            atname = normalize_hydrogen_name(atom.get_name())
                            atom_dict[atomidx] = atom
                            coord = atom.get_coord()
                            x = coord[0]
                            y = coord[1]
                            z = coord[2]
                            pl_set_atom(PROTEIN, atomidx, x, y, z, residx, chainidx, atname)
                            atomidx = atomidx + 1
                    atomidx1 = atomidx - 1

                    firstatom = firstatom + [atomidx0]
                    lastatom = lastatom + [atomidx1]

                    pl_set_residue(PROTEIN, residx, firstatom[residx], lastatom[residx], chainidx, aaname)
                    residx = residx + 1
                    
            lastres = lastres + [residx - 1]
            pl_set_chain(PROTEIN, chainidx, firstres[chainidx], lastres[chainidx])
            chainidx = chainidx + 1
        print "DONE."

        # Building protein connectivity based on the atomic data just copied.
        print "BUILDING MOLECULE CONNECTIVITY..."
        pl_build_connectivity(PROTEIN, use_sc)
        print "DONE."

        # Running simulation with curent atomic information and simulation driver.
        print "RUNNING SIMULATION..."
        pl_run_simulation(cfg)
        print "SIMULATION COMPLETED."

        print "COPYING DATA BACK TO BIOPYTHON...",
        # Copy coordinates back to the stucture object.
        for atomkey in atom_dict.keys():
            newcoord = pl_get_atom_coords(atomkey)
            atom = atom_dict[atomkey]
            coord = atom.get_coord()
            coord[0] = newcoord[0]
            coord[1] = newcoord[1]
            coord[2] = newcoord[2]
        print "DONE."
      
        #print atom_dict[50].get_parent(),atom_dict[50],atom_dict[50].get_coord()
        #for chain in model: 
        #    for residue in chain:
        #        if residue.get_id()[1]==11:
        #            for atom in residue:
        #                print atom,atom.get_coord()

        #for chain in self.proteins[0]: 
        #    for residue in chain:
        #        if residue.get_id()[1]==11:
        #            for atom in residue:
        #                print atom,atom.get_coord()

        # Deleting the data arrays in C.
        print "RELEASING MEMORY USED IN SIMULATION...",
        pl_free_memory()
        print "DONE."
