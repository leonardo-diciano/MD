################################################################
# Script to calculate Mulliken Charges from a HF/STO-3G calculation in PSI4
# charges are printed to the terminal
# if used with smiles, the corresponding xyzfile will be created
# limited to elements until Ar; due to manual definition of basis set /electron numbers
# 
# Author: Lila Zapp (2025)
#
################################################################
# use this script with smiles, e.g.:
# python3 get_mulliken_charges.py -smi "CCC"
# or with xyz file:
# python3 get_mulliken_charges.py -xyz propane.xyz
################################################################

import sys
import os
import numpy as np
import psi4
from rdkit import Chem
from rdkit.Chem import AllChem

np.set_printoptions(precision=8, linewidth=200)


if sys.argv[1] == '-smi':
    #make moldoc from smiles
    mol=Chem.MolFromSmiles(sys.argv[2])
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    moldoc=Chem.MolToXYZBlock(mol)

    # generate xyz file
    xyzfile = sys.argv[2]+".xyz"
    if os.path.exists(xyzfile):
        print("%s will be replaced"%xyzfile)
        os.remove(xyzfile)
    with open(xyzfile,'w') as file:
        file.write(moldoc)
    file.close()
    print("wrote smiles to %s file"%xyzfile)
elif sys.argv[1] == '-xyz':
    #make moldoc from xyz file
    with open(str(sys.argv[2])) as file:
        moldoc = ""
        for line in file:
            moldoc += line
else:
    print("Use either the -xyz or -smi flag")
    exit


#print(moldoc)


first_row = ["H","He"]
second_row_s = ["Li","Be"]
second_row_p = ["B","C","N","O","F","Ne"]
third_row_s = ["Na","Mg"]
third_row_p = ["Al","Si","P","S","Cl","Ar"]

element_names = []
for element in first_row:
    element_names.append(element)
for element in second_row_s:
    element_names.append(element)
for element in second_row_p:
    element_names.append(element)
for element in third_row_s:
    element_names.append(element)
for element in third_row_p:
    element_names.append(element)


# assign number of electrons
n_electrons = {}
atomic_no = 1
for element in element_names:
    n_electrons[element] = atomic_no
    atomic_no += 1

#print(n_electrons)


# assign number of basis functions per atom in STO-3G
nbf_sto3g = {}
for element in element_names:
    if element in first_row:
        nbf_sto3g[element] = 1
    elif element in second_row_s:
        nbf_sto3g[element] = 2
    elif element in second_row_p:
        nbf_sto3g[element] = 5
    elif element in third_row_s:
        nbf_sto3g[element] = 6
    elif element in third_row_p:
        nbf_sto3g[element] = 9

#print(nbf_sto3g)



#### The Mulliken Charges are also part of the molden file, so let us calculate them:
def get_density_mat(C,n_occ):
    nbf = C.shape[0]
    D = np.zeros((nbf,nbf))
    for mu in range(nbf):
        for nu in range(nbf):
            for i in range(n_occ):
                D[mu,nu] += 2*C[mu,i]*C[nu,i]
    return D

def get_Mulliken_charge(C,S,n_occ,nbf_on_A, nel_on_A):
    #calculate density matrix
    D = get_density_mat(C=C,n_occ=n_occ)

    #calculate population matrix (elementwise multiplication of D and S)
    #population = np.zeros((nbf,nbf))

    #Pmunu = Dmunu * Smunu
    population = D * S

    #get gross orbital product for orbital mu (GOP_mu) by summing over nu
    nbf = population.shape[0]
    totalGOP=0 #=N_elec
    GOP_mus = []
    for mu in range(nbf):
        GOP_mu = 0
        for nu in range(nbf):
            GOP_mu += population[mu,nu]
        totalGOP += GOP_mu
        GOP_mus.append(GOP_mu)

    GOP_mus = np.asarray(GOP_mus)

    index = 0
    GOP_A = np.zeros(len(nbf_on_A))
    for atom in range(len(nbf_on_A)):
        #print(index,nbf_on_A[atom]+nbf_on_A[atom])
        for mu in range(index,index+nbf_on_A[atom]):
            GOP_A[atom] += GOP_mus[mu]
        index += nbf_on_A[atom]

    #print('GOP_mus =           ',GOP_mus)
    #print('GOP_A =             ',GOP_A)
    #print('total GOP = N_elec =',totalGOP)
    return nel_on_A - GOP_A


psi4.core.set_output_file('PSI4_output.dat',False)

mol=psi4.geometry(moldoc)
mol.update_geometry()


# DEFINE nbf_on_A and nel_on_A
#assign for xyz
atom_list = []
for line in moldoc.split("\n"):
    for element in element_names:
        if line.startswith(element):
            atom_list.append(element)
#print(atom_list)

nbf_on_A = []
for atom in atom_list:
    nbf_on_A.append(nbf_sto3g[atom])
nbf_on_A = np.asarray(nbf_on_A)

nel_on_A = []
for atom in atom_list:
    nel_on_A.append(n_electrons[atom])
nel_on_A = np.asarray(nel_on_A)


#e.g. for propane:
#nbf_on_A = np.array([5,5,5,1,1,1,1,1,1,1,1])
#nel_on_A = np.array([6,6,6,1,1,1,1,1,1,1,1])
#print(nbf_on_A,nel_on_A)

#specification of a basis set and additional options
psi4.set_options({'basis':        'sto-3g',
                  'scf_type':         'pk',
                  'e_convergence':    1e-8,
                  'd_convergence':    1e-8})

#Hartree-Fock calculation for the water molecule
energy,wfn=psi4.energy('hf',molecule=mol,return_wfn=True)

mints = psi4.core.MintsHelper(wfn.basisset())

#print("HF energy of the water molecule:",energy)
CMO = np.asarray(wfn.Ca())
S = np.asarray(mints.ao_overlap())
n_occ = wfn.nalpha()


if n_occ != sum(nel_on_A)/2:
    print("WARNING: check no of electrons, sth doesnt add up")
if S.shape[0] != sum(nbf_on_A):
    print("WARNING: check no of basis functions, sth doesnt add up")


print("\nMulliken charges: ",get_Mulliken_charge(CMO,S,n_occ,nbf_on_A, nel_on_A))
