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

# or if you want to print the charges to an extra file: add the "-f" flag in the end for "file"
# python3 get_mulliken_charges.py -smi "CCC" -f
################################################################


import numpy as np
import psi4
from rdkit import Chem
from rdkit.Chem import AllChem

np.set_printoptions(precision=8, linewidth=200)

class Mulliken_Charges():
    def __init__(self,xyzblock):
        self.moldoc = "\n".join(xyzblock.splitlines()[2:])
        self.initialize()
        return 
    
    def initialize(self):
        first_row = ["H","He"]
        second_row_s = ["Li","Be"]
        second_row_p = ["B","C","N","O","F","Ne"]
        third_row_s = ["Na","Mg"]
        third_row_p = ["Al","Si","P","S","Cl","Ar"]

        self.element_names = []
        for element in first_row:
            self.element_names.append(element)
        for element in second_row_s:
            self.element_names.append(element)
        for element in second_row_p:
            self.element_names.append(element)
        for element in third_row_s:
            self.element_names.append(element)
        for element in third_row_p:
            self.element_names.append(element)


        # assign number of electrons
        self.n_electrons = {}
        atomic_no = 1
        for element in self.element_names:
            self.n_electrons[element] = atomic_no
            atomic_no += 1

        #print(self.n_electrons)


        # assign number of basis functions per atom in STO-3G
        self.nbf_sto3g = {}
        for element in self.element_names:
            if element in first_row:
                self.nbf_sto3g[element] = 1
            elif element in second_row_s:
                self.nbf_sto3g[element] = 2
            elif element in second_row_p:
                self.nbf_sto3g[element] = 5
            elif element in third_row_s:
                self.nbf_sto3g[element] = 6
            elif element in third_row_p:
                self.nbf_sto3g[element] = 9

    #### The Mulliken Charges are also part of the molden file, so let us calculate them:
    def get_density_mat(self,C,n_occ):
        nbf = C.shape[0]
        D = np.zeros((nbf,nbf))
        for mu in range(nbf):
            for nu in range(nbf):
                for i in range(n_occ):
                    D[mu,nu] += 2*C[mu,i]*C[nu,i]
        return D

    def get_Mulliken_charge(self,C,S,n_occ,nbf_on_A, nel_on_A):
        #calculate density matrix
        D = self.get_density_mat(C=C,n_occ=n_occ)

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
    
    def main(self):
        psi4.core.set_output_file('PSI4_output.dat',False)

        mol=psi4.geometry(f"""
                  symmetry c1
                  {self.moldoc}
                  """)
        mol.update_geometry()

        # DEFINE nbf_on_A and nel_on_A
        #assign for xyz
        atom_list = []
        for line in self.moldoc.split("\n"):
            for element in self.element_names:
                if line.startswith(element):
                    atom_list.append(element)
        #print(atom_list)

        nbf_on_A = []
        for atom in atom_list:
            nbf_on_A.append(self.nbf_sto3g[atom])
        nbf_on_A = np.asarray(nbf_on_A)

        nel_on_A = []
        for atom in atom_list:
            nel_on_A.append(self.n_electrons[atom])
        nel_on_A = np.asarray(nel_on_A)


        #e.g. for propane:
        #nbf_on_A = np.array([5,5,5,1,1,1,1,1,1,1,1])
        #nel_on_A = np.array([6,6,6,1,1,1,1,1,1,1,1])
        #print(nbf_on_A,nel_on_A)

        #specification of a basis set and additional options
        psi4.set_options({'basis':        'sto-3g',
                        'scf_type':         'pk',
                        'e_convergence':    1e-8,
                        'd_convergence':    1e-8,})

        #Hartree-Fock calculation for the water molecule
        energy,wfn=psi4.energy('hf',molecule=mol,return_wfn=True)

        mints = psi4.core.MintsHelper(wfn.basisset())

        #print("HF energy of the water molecule:",energy)
        CMO = np.asarray(wfn.Ca())
        S = np.asarray(mints.ao_overlap())
        n_occ = wfn.nalpha()


        if n_occ != sum(nel_on_A)/2:
            raise ValueError("WARNING: check no of electrons, sth doesnt add up")
        if S.shape[0] != sum(nbf_on_A):
            raise ValueError("WARNING: check no of basis functions, sth doesnt add up")


        charges = self.get_Mulliken_charge(CMO,S,n_occ,nbf_on_A, nel_on_A)
        return charges

# write to file
#if len(sys.argv) > 3:
#    if sys.argv[3] == '-f':
#        with open(sys.argv[2]+"_charges.txt","w") as file:        
#            file.write("[Charges]\n")
#            natom=0
#            for charge in charges:
#                natom +=1            
#                file.write("%3i     %16.12f\n"%(natom,charge))
#        file.close()
#        print("wrote Mulliken charges to %s"%(sys.argv[2]+"_charges.txt"))
#    else:
#        print("\nMulliken charges: ",charges)




