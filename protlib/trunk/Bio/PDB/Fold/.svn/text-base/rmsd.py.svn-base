#!/usr/bin/env python
from Bio.PDB import *
from Bio import pairwise2
import sys

__doc__="This module contains the functionality to calculate the rmsd between two three-dimensional structures"



def align_structures(protein1,protein2,model1=0,model2=0):
    """This function takes two structures as pdb parser objects and returns an alignment of the two structures"""
    seq1list=[x.get_sequence() for x in PPBuilder().build_peptides(protein1)]
    seq2list=[x.get_sequence() for x in PPBuilder().build_peptides(protein2)]

    seq1=seq1list[model1]
    seq2=seq2list[model2]
    return pairwise2.align.globalxx(str(seq1),str(seq2))

def calculate_rmsd(structure1,structure2,model1=0,model2=0,uselist=['CA']):
    """This is the main function for calculating the rmsd between two structures"""
    if hasattr(structure1,'get_chains'): #already a pdb parser object
        protein1=structure1 
    else: 
        protein1=PDBParser().get_structure("protein1",structure1)
    if hasattr(structure2,'get_chains'): #already a pdb parser object
        protein2=structure2 
    else: 
        protein2=PDBParser().get_structure("protein2",structure2)

    alignment=align_structures(protein1,protein2,model1,model2)
    if len(alignment)>2: alignment=alignment[2]
    else: alignment=alignment[0]

    residues_1=[x for x in protein1.get_residues() ]
    residues_2=[x for x in protein2.get_residues() ]
    fixed_list=[]
    moving_list=[]
    count=0
    for resnum in range(len(alignment[0])):
        res=alignment[0][resnum]
        res2=alignment[1][resnum]
        if res in Polypeptide.aa1 and res2 in Polypeptide.aa1:
            for atom in residues_1[count]:
                if atom.get_name() in uselist:
                    fixed_list.append(atom)
            count=count+1
    count=0
    for resnum in range(len(alignment[1])):
        res1=alignment[0][resnum]
        res=alignment[1][resnum]
        
        if res in Polypeptide.aa1 and res1 in Polypeptide.aa1:
            for atom in residues_2[count]:
                if atom.get_name() in uselist:
                    moving_list.append(atom)
            count=count+1
       
    sup=Superimposer()
    sup.set_atoms(fixed_list,moving_list)
    sup.apply(protein2.get_atoms())

    return sup.rms,protein2

if __name__=="__main__":
    structure1=sys.argv[1]
    structure2=sys.argv[2]
    print calculate_rmsd(structure1,structure2)[0]
