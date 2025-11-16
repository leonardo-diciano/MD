"""Atom Typer and Topology File Generator
Starting from SMILES (python get_topology.py -smi 'SMILES' output_filename.top) or
XYZ (python get_topology.py -xyz coord_file.xyz output_filename.top), it assigns the atom types
according to definitions in force_field_data.py. Then, it proceeds to scan all the molecule bonds,
angles, proper and improper dihedral to assign the corresponding parameters to each and store them in lists. 
It also assigns parameteres for Lennard-Jones potential and it calculates the RESP charges for each
atom at HF/STO-3G level using psi4 and its resp plugin. Finally, the whole set of parameters needed for the 
chosen molecule are printed in a *.top file, ordered by type of interaction.

The atom types are a reduced version of AMBER14sb_parmbsc1, as well as all parameters and are taken from:
https://github.com/intbio/gromacs_ff/tree/master/amber14sb_parmbsc1_opc_lmi.ff .

In the output file, atoms' index count starts at 1, instead of 0, to improve compatibility with Fortran scripts.
Example: python get_topology.py -smi 'c1ccccc1' 'benzene.top'

Author: Leonardo Di Ciano (2025)"""

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from force_field_data import atom_types, bond_params, angle_params, improper_dihedral_params, dihedral_params, LJ_params
import sys
import itertools
import warnings
import psi4
import resp
import os
import glob

#Avoid pandas FutureWarnings
warnings.filterwarnings('ignore',category=FutureWarning)


#Import atom types and MW from dict
data=pd.DataFrame.from_dict(atom_types,orient='index',columns=["AtomType","MW"]) 
atom_types=data["AtomType"].to_list() 
atom_types=[x.strip() for x in atom_types]

#SMILES or XYZ file run mode
if sys.argv[1] == '-smi':
    mol=Chem.MolFromSmiles(sys.argv[2])
    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
elif sys.argv[2] == '-xyz':
    mol=Chem.MolFromXYZFile(sys.argv[2])
    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
else:
    print("Use either the -xyz or -smi flag")
    exit

"""Atom Types Definition
The assignment is generally based on combinations of type of neighbor atoms, 
atoms hybridization and aromatic nature."""
atom_assigned_types={} #dict with index and atom types
list_atom_types_to_print=[] #list with tuples of atom idx (+1), symbol, type and MW
for atom_i in mol.GetAtoms():
    if atom_i.GetSymbol() == "C":
        symbols=[]
        for i in atom_i.GetNeighbors():
            symbols.append(i.GetSymbol())

        if "O" in symbols and atom_i.GetHybridization() == "SP2":
            atom_assigned_types[atom_i.GetIdx()]="C"
            list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"C",
                                                float(data["MW"][data["AtomType"]=="C"])))
        elif atom_i.GetIsAromatic() == True or atom_i.GetHybridization() == "SP2":
            atom_assigned_types[atom_i.GetIdx()]="CA"
            list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"CA",
                                                float(data["MW"][data["AtomType"]=="CA"])))
        else:
            atom_assigned_types[atom_i.GetIdx()]="CT"
            list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"CT",
                                                float(data["MW"][data["AtomType"]=="CT"])))
        continue
    elif atom_i.GetSymbol() == "H":
        for nghb in atom_i.GetNeighbors():
            if nghb.GetSymbol() == "C" and nghb.GetIsAromatic() == True:
                atom_assigned_types[atom_i.GetIdx()]="HA"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"HA",
                                                float(data["MW"][data["AtomType"]=="HA"])))
            elif nghb.GetSymbol() == "C":
                atom_assigned_types[atom_i.GetIdx()]="HC"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"HC",
                                                float(data["MW"][data["AtomType"]=="HC"])))
            elif nghb.GetSymbol() == "O":
                atom_assigned_types[atom_i.GetIdx()]="HO"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"HO",
                                                float(data["MW"][data["AtomType"]=="HO"])))
            elif nghb.GetSymbol() == "N":
                atom_assigned_types[atom_i.GetIdx()]="H"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"H",
                                                float(data["MW"][data["AtomType"]=="H"])))    
            elif nghb.GetSymbol() == "S":
                atom_assigned_types[atom_i.GetIdx()]="HS"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"HS",
                                                float(data["MW"][data["AtomType"]=="HS"])))
        continue 
    
    elif atom_i.GetSymbol() == "N":
        if len(atom_i.GetNeighbors()) == 2:
            for nghb in atom_i.GetNeighbors():
                if nghb.GetSymbol() == "H":
                    continue
                elif nghb.GetSymbol() == "C" and nghb.GetHybridization() == "SP2" and atom_i.GetHybridization() == "SP2" :
                    atom_assigned_types[atom_i.GetIdx()]="N"
                    list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"N",
                                                    float(data["MW"][data["AtomType"]=="N"]))) 
                else :
                    atom_assigned_types[atom_i.GetIdx()]="N*"
                    list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"N*",
                                                    float(data["MW"][data["AtomType"]=="N*"]))) 
        elif len(atom_i.GetNeighbors()) == 3:
            for nghb in atom_i.GetNeighbors():
                if nghb.GetSymbol() == "H":
                    continue
                elif nghb.GetSymbol() == "C" and atom_i.GetHybridization() == "SP2" :
                    atom_assigned_types[atom_i.GetIdx()]="N2"
                    list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"N2",
                                                    float(data["MW"][data["AtomType"]=="N2"])))
                else:
                    atom_assigned_types[atom_i.GetIdx()]="N*"
                    list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"N*",
                                                    float(data["MW"][data["AtomType"]=="N*"])))
        elif len(atom_i.GetNeighbors()) == 4:
            atom_assigned_types[atom_i.GetIdx()]="N3"
            list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"N3",
                                            float(data["MW"][data["AtomType"]=="N3"])))
        else:
            atom_assigned_types[atom_i.GetIdx()]="N*"
            list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"N*",
                                            float(data["MW"][data["AtomType"]=="N*"])))
        continue
    elif atom_i.GetSymbol() == "O":
        if len(atom_i.GetNeighbors()) == 1:
            for nghb in atom_i.GetNeighbors():
                if nghb.GetSymbol() == "H":
                    atom_assigned_types[atom_i.GetIdx()]="OH"
                    list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"OH",
                                                    float(data["MW"][data["AtomType"]=="OH"]))) 
                elif nghb.GetSymbol() == "C" and atom_i.GetFormalCharge() == -1:
                    atom_assigned_types[atom_i.GetIdx()]="O2"
                    list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"O2",
                                                    float(data["MW"][data["AtomType"]=="O2"]))) 
                elif nghb.GetSymbol() == "C" and nghb.GetHybridization() == "SP2":
                    atom_assigned_types[atom_i.GetIdx()]="O"
                    list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"O",
                                                    float(data["MW"][data["AtomType"]=="O"]))) 
                else :
                    atom_assigned_types[atom_i.GetIdx()]="OP"
                    list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"OP",
                                                    float(data["MW"][data["AtomType"]=="OP"]))) 
        elif len(atom_i.GetNeighbors()) == 2:
            nghb1=atom_i.GetNeighbors()[0]
            nghb2=atom_i.GetNeighbors()[1]
            if nghb1.GetSymbol() == "H" or nghb2.GetSymbol() == "H":
                atom_assigned_types[atom_i.GetIdx()]="OH"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"OH",
                                                float(data["MW"][data["AtomType"]=="OH"])))
            elif nghb1.GetSymbol() == "C" or nghb2.GetSymbol() == "C":
                atom_assigned_types[atom_i.GetIdx()]="OS"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"OS",
                                                float(data["MW"][data["AtomType"]=="OS"])))
            else :
                atom_assigned_types[atom_i.GetIdx()]="OP"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"OP",
                                                float(data["MW"][data["AtomType"]=="OP"])))
        continue 
    elif atom_i.GetSymbol() == "S":
        for nghb in atom_i.GetNeighbors():
            if nghb.GetSymbol() == "H":
                atom_assigned_types[atom_i.GetIdx()]="SH"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"SH",
                                                float(data["MW"][data["AtomType"]=="SH"])))
            else:
                atom_assigned_types[atom_i.GetIdx()]="S"
                list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),"S",
                                                float(data["MW"][data["AtomType"]=="S"])))
        continue
    else:
        atom_assigned_types[atom_i.GetIdx()]=atom_i.GetSymbol()
        list_atom_types_to_print.append((atom_i.GetIdx()+1,atom_i.GetSymbol(),atom_i.GetSymbol(),
                                         float(data["MW"][data["AtomType"]==atom_i.GetSymbol()])))
        continue

print(f"Atom Types summary\n",list_atom_types_to_print)

"""Bonds assignments
Assigned simply by iterating over the bonds from rdkit and getting corresponding params
from the force_field_data library"""
list_bond_params_to_print=[]
for bond_i in mol.GetBonds():
    a1=bond_i.GetBeginAtomIdx()
    a2=bond_i.GetEndAtomIdx()
    par_12=bond_params[f"{atom_assigned_types[a1]}-{atom_assigned_types[a2]}"]
    list_bond_params_to_print.append((a1+1,a2+1,par_12[0],par_12[1]))

print(f"Bonding summary\n",list_bond_params_to_print)

"""Angles assignments
Assigned by iterating over the possible combinations of a central atom and two *different* neighbors,
itertools.combination avoids repetitions (ie. for three neighbors, pairs contain [['A','B'],['A','C'],['B','C']] 
but not ['B','A'] too). The parameters are then collected from the force_field_data library"""
list_angle_params_to_print=[]
for atom_i in mol.GetAtoms():
    neighbors=atom_i.GetNeighbors()
    if len(atom_i.GetNeighbors()) >= 2:
        symbols=[]
        for i in atom_i.GetNeighbors():
            symbols.append(atom_assigned_types[i.GetIdx()])
        pairs = list(itertools.combinations(symbols, 2))
        for pair in pairs:
           aa=f"{pair[0]}-{atom_assigned_types[atom_i.GetIdx()]}-{pair[1]}"
           par_aa=angle_params[aa]
           list_angle_params_to_print.append((aa,pair[0]+1,atom_i.GetIdx()+1,pair[1]+1,par_aa[0],par_aa[1])) 
    else:
        continue

print(f"Angles summary\n",list_angle_params_to_print)

"""Improper Diehedrals assignments
Assigned by iterating over the possible combinations of a central atom and three *different* neighbors,
itertools.permutations allows repetitions, in order to get a correctly ordered set to grep from the 
force_field_data library"""
list_impro_dihedrals_params_to_print=[]
for atom_i in mol.GetAtoms():
    neighbors=atom_i.GetNeighbors()
    if len(atom_i.GetNeighbors()) == 3:
        symbols=[]
        for i in atom_i.GetNeighbors():
            symbols.append(atom_assigned_types[i.GetIdx()])
        tris = list(itertools.permutations(symbols, 3))
        for tri in tris:
            aa=f"{atom_assigned_types[atom_i.GetIdx()]}-{tri[0]}-{tri[1]}-{tri[2]}"
            if aa in improper_dihedral_params.keys():
                par_aa=improper_dihedral_params[aa]
                list_impro_dihedrals_params_to_print.append((aa,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_aa[0],par_aa[1],par_aa[2])) 
            else:
                continue
    else:
        continue
print(f"Improper dihedrals summary\n",list_impro_dihedrals_params_to_print)

"""Proper Diehedrals assignments
Assigned by iterating over the neigbors(nghb) of the first (i) and second (j) central atoms, 
until obtaining a i_nghb-i-j-j_nghb or j_nghb-j-i-i_nghb combination that matches one of the 
force_field_data library entries."""
list_dihedrals_params_to_print=[]
for atom_i in mol.GetAtoms():
    neighbors=atom_i.GetNeighbors()
    if len(neighbors) < 2:
        continue
    else:
        for atom_j in neighbors:
            if len(atom_j.GetNeighbors()) <= 1:
                continue
            else:
                for j_nghb in atom_j.GetNeighbors():
                    if j_nghb.GetIdx() != atom_i.GetIdx():
                        for i_nghb in neighbors:
                            if i_nghb.GetIdx() != atom_j.GetIdx():
                                aa=f"{atom_assigned_types[i_nghb.GetIdx()]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[atom_j.GetIdx()]}-{atom_assigned_types[j_nghb.GetIdx()]}"
                                bb=f"{atom_assigned_types[j_nghb.GetIdx()]}-{atom_assigned_types[atom_j.GetIdx()]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[i_nghb.GetIdx()]}"
                                if aa in dihedral_params.keys():
                                    par_aa=dihedral_params[aa]
                                    list_dihedrals_params_to_print.append((aa,i_nghb.GetIdx()+1,atom_i.GetIdx()+1,atom_j.GetIdx()+1,j_nghb.GetIdx()+1,par_aa[0],par_aa[1],par_aa[2]))
                                elif bb in dihedral_params.keys():
                                    par_bb=dihedral_params[bb]
                                    list_dihedrals_params_to_print.append((bb,j_nghb.GetIdx()+1,atom_j.GetIdx()+1,atom_i.GetIdx()+1,i_nghb.GetIdx()+1,par_bb[0],par_bb[1],par_bb[2]))
                                else:
                                    continue
                            else:
                                continue
                        else:
                            continue
                    else:
                        continue

print(f"Proper dihedrals summary\n",list_dihedrals_params_to_print)

"""Lennard-Jones assignments
The parameters are simply assigned for each atom type."""
list_LJ_params_to_print=[]
for atom_i in mol.GetAtoms():
    idx=atom_i.GetIdx()
    type_i=atom_assigned_types[idx]
    par_lj=LJ_params[type_i]
    list_LJ_params_to_print.append((idx+1,par_lj[0],par_lj[1]))

print(f"LJ params summary\n",list_LJ_params_to_print)

"""RESP charges assignments
The charges for Coulomb electrostatic interaction are not available in the library, but are calculated as 
RESP charges. In this case, it uses the resp plugin of psi4 and there charges are obtained with two iterations,
following https://github.com/cdsgroup/resp/blob/master/examples/example1.py """
list_resp_charges_to_print=[]
print("Calculating RESP Charges with Psi4's RESP plugin")
psi_mol = psi4.geometry(Chem.MolToXYZBlock(mol))
psi_mol.update_geometry()
options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
           }
# Call for first stage fit
charges1 = resp.resp([psi_mol], options)
# Change the value of the RESP parameter A
options['RESP_A'] = 0.001
# Add constraint for atoms fixed in second stage fit
constraint_charge = []
for i in range(4, 8):
    constraint_charge.append([charges1[1][i], [i+1]])
options['constraint_charge'] = constraint_charge
options['constraint_group'] = [[2, 3, 4]]
options['grid'] = ['1_%s_grid.dat' %psi_mol.name()]
options['esp'] = ['1_%s_grid_esp.dat' %psi_mol.name()]
# Call for second stage fit
charges2 = resp.resp([psi_mol], options)
# Get RESP charges
for i in range(len(charges2[1])):
    list_resp_charges_to_print.append((i+1,charges2[1][i]))

for i in glob.glob("*.dat"):
    os.remove(i)
os.remove("results.out")

print('RESP Charges Summary:',list_resp_charges_to_print)

"""The parameters are assembled with all the stored components and
written in a *.top file"""
if sys.argv[3]:
    file=sys.argv[3]
else:
    file="topology.top"
with open(file,"w+") as f:
    f.write(f"[AtomTypes] {len(list_atom_types_to_print)}\n")
    types_df=pd.DataFrame(list_atom_types_to_print,columns=["Idx","Symb","Type","MW"],index=None)
    types_df.to_csv(f,sep="\t",header=False,index=False)
    f.write(f"\n[Bonds] {len(list_bond_params_to_print)}\n")
    bonds_df=pd.DataFrame(list_bond_params_to_print,columns=["Idx1","Idx2","b0","kb"],index=None)
    bonds_df.to_csv(f,sep="\t",header=False,index=False)
    f.write(f"\n[Angles] {len(list_angle_params_to_print)}\n")
    angles_df=pd.DataFrame(list_angle_params_to_print,columns=["Types","Idx1","Idx2","Idx3","th0","cth"],index=None)
    angles_df.to_csv(f,sep="\t",header=False,index=False)
    f.write(f"\n[ImproperDihedrals] {len(list_impro_dihedrals_params_to_print)}\n")
    impdie_df=pd.DataFrame(list_impro_dihedrals_params_to_print,columns=["Types","Idx1","Idx2","Idx3","Idx4","phase","kd","pn"],index=None)
    impdie_df.to_csv(f,sep="\t",header=False,index=False)
    f.write(f"\n[Dihedrals] {len(list_dihedrals_params_to_print)}\n")
    die_df=pd.DataFrame(list_dihedrals_params_to_print,columns=["Types","Idx1","Idx2","Idx3","Idx4","phase","kd","pn"],index=None)
    die_df.to_csv(f,sep="\t",header=False,index=False)
    f.write(f"\n[LJ]\n")
    lj_df=pd.DataFrame(list_LJ_params_to_print,columns=["Idx","sigma","epsilon"],index=None)
    lj_df.to_csv(f,sep="\t",header=False,index=False)
    f.write(f"\n[Charges]\n")
    chg_df=pd.DataFrame(list_resp_charges_to_print,columns=["Idx","charge"],index=None)
    chg_df.to_csv(f,sep="\t",header=False,index=False)
    f.write(f" ")
print(f"Topology file written in {file}")
