from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from force_field_data import atom_types, bond_params, angle_params, improper_dihedral_params, dihedral_params, LJ_params
import sys
import itertools

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


#Atom Types Definition
atom_assigned_types={} #dict with index and atom types
list_atom_types_to_print=[] #list with tuples of atom idx, symbol, type and MW
for atom_i in mol.GetAtoms():
    if atom_i.GetSymbol() == "C":
        symbols=[]
        for i in atom_i.GetNeighbors():
            symbols.append(i.GetSymbol())

        if "O" in symbols and atom_i.GetHybridization() == "SP2":
            atom_assigned_types[atom_i.GetIdx()]="C"
            list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"C",
                                                float(data["MW"][data["AtomType"]=="C"])))
        elif atom_i.GetIsAromatic() == True or atom_i.GetHybridization() == "SP2":
            atom_assigned_types[atom_i.GetIdx()]="CA"
            list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"CA",
                                                float(data["MW"][data["AtomType"]=="CA"])))
        else:
            atom_assigned_types[atom_i.GetIdx()]="CT"
            list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"CT",
                                                float(data["MW"][data["AtomType"]=="CT"])))
        continue
    elif atom_i.GetSymbol() == "H":
        for nghb in atom_i.GetNeighbors():
            if nghb.GetSymbol() == "C" and nghb.GetIsAromatic() == True:
                atom_assigned_types[atom_i.GetIdx()]="HA"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"HA",
                                                float(data["MW"][data["AtomType"]=="HA"])))
            elif nghb.GetSymbol() == "C":
                atom_assigned_types[atom_i.GetIdx()]="HC"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"HC",
                                                float(data["MW"][data["AtomType"]=="HC"])))
            elif nghb.GetSymbol() == "O":
                atom_assigned_types[atom_i.GetIdx()]="HO"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"HO",
                                                float(data["MW"][data["AtomType"]=="HO"])))
            elif nghb.GetSymbol() == "N":
                atom_assigned_types[atom_i.GetIdx()]="H"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"H",
                                                float(data["MW"][data["AtomType"]=="H"])))    
            elif nghb.GetSymbol() == "S":
                atom_assigned_types[atom_i.GetIdx()]="HS"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"HS",
                                                float(data["MW"][data["AtomType"]=="HS"])))
        continue 
    
    elif atom_i.GetSymbol() == "N":
        if len(atom_i.GetNeighbors()) == 2:
            for nghb in atom_i.GetNeighbors():
                if nghb.GetSymbol() == "H":
                    continue
                elif nghb.GetSymbol() == "C" and nghb.GetHybridization() == "SP2" and atom_i.GetHybridization() == "SP2" :
                    atom_assigned_types[atom_i.GetIdx()]="N"
                    list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"N",
                                                    float(data["MW"][data["AtomType"]=="N"]))) 
                else :
                    atom_assigned_types[atom_i.GetIdx()]="N*"
                    list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"N*",
                                                    float(data["MW"][data["AtomType"]=="N*"]))) 
        elif len(atom_i.GetNeighbors()) == 3:
            for nghb in atom_i.GetNeighbors():
                if nghb.GetSymbol() == "H":
                    continue
                elif nghb.GetSymbol() == "C" and atom_i.GetHybridization() == "SP2" :
                    atom_assigned_types[atom_i.GetIdx()]="N2"
                    list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"N2",
                                                    float(data["MW"][data["AtomType"]=="N2"])))
                else:
                    atom_assigned_types[atom_i.GetIdx()]="N*"
                    list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"N*",
                                                    float(data["MW"][data["AtomType"]=="N*"])))
        elif len(atom_i.GetNeighbors()) == 4:
            atom_assigned_types[atom_i.GetIdx()]="N3"
            list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"N3",
                                            float(data["MW"][data["AtomType"]=="N3"])))
        else:
            atom_assigned_types[atom_i.GetIdx()]="N*"
            list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"N*",
                                            float(data["MW"][data["AtomType"]=="N*"])))
        continue
    elif atom_i.GetSymbol() == "O":
        if len(atom_i.GetNeighbors()) == 1:
            for nghb in atom_i.GetNeighbors():
                if nghb.GetSymbol() == "H":
                    atom_assigned_types[atom_i.GetIdx()]="OH"
                    list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"OH",
                                                    float(data["MW"][data["AtomType"]=="OH"]))) 
                elif nghb.GetSymbol() == "C" and atom_i.GetFormalCharge() == -1:
                    atom_assigned_types[atom_i.GetIdx()]="O2"
                    list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"O2",
                                                    float(data["MW"][data["AtomType"]=="O2"]))) 
                elif nghb.GetSymbol() == "C" and nghb.GetHybridization() == "SP2":
                    atom_assigned_types[atom_i.GetIdx()]="O"
                    list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"O",
                                                    float(data["MW"][data["AtomType"]=="O"]))) 
                else :
                    atom_assigned_types[atom_i.GetIdx()]="OP"
                    list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"OP",
                                                    float(data["MW"][data["AtomType"]=="OP"]))) 
        elif len(atom_i.GetNeighbors()) == 2:
            nghb1=atom_i.GetNeighbors()[0]
            nghb2=atom_i.GetNeighbors()[1]
            if nghb1.GetSymbol() == "H" or nghb2.GetSymbol() == "H":
                atom_assigned_types[atom_i.GetIdx()]="OH"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"OH",
                                                float(data["MW"][data["AtomType"]=="OH"])))
            elif nghb1.GetSymbol() == "C" or nghb2.GetSymbol() == "C":
                atom_assigned_types[atom_i.GetIdx()]="OS"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"OS",
                                                float(data["MW"][data["AtomType"]=="OS"])))
            else :
                atom_assigned_types[atom_i.GetIdx()]="OP"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"OP",
                                                float(data["MW"][data["AtomType"]=="OP"])))
        continue 
    elif atom_i.GetSymbol() == "S":
        for nghb in atom_i.GetNeighbors():
            if nghb.GetSymbol() == "H":
                atom_assigned_types[atom_i.GetIdx()]="SH"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"SH",
                                                float(data["MW"][data["AtomType"]=="SH"])))
            else:
                atom_assigned_types[atom_i.GetIdx()]="S"
                list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),"S",
                                                float(data["MW"][data["AtomType"]=="S"])))
        continue
    else:
        atom_assigned_types[atom_i.GetIdx()]=atom_i.GetSymbol()
        list_atom_types_to_print.append((atom_i.GetIdx(),atom_i.GetSymbol(),atom_i.GetSymbol(),
                                         float(data["MW"][data["AtomType"]==atom_i.GetSymbol()])))
        continue

print(f"Atom Types summary\n",list_atom_types_to_print)

#Bonds assignments

list_bond_params_to_print=[]
for bond_i in mol.GetBonds():
    a1=bond_i.GetBeginAtomIdx()
    a2=bond_i.GetEndAtomIdx()
    par_12=bond_params[f"{atom_assigned_types[a1]}-{atom_assigned_types[a2]}"]
    list_bond_params_to_print.append((a1,a2,par_12[0],par_12[1]))

print(f"Bonding summary\n",list_bond_params_to_print)

#Angles assignments
list_angle_params_to_print=[]
for atom_i in mol.GetAtoms():
    neighbors=atom_i.GetNeighbors()
    if atom_i.GetNeighbors() >= 2:
        symbols=[]
        for i in atom_i.GetNeighbors():
            symbols.append(atom_assigned_types[i.GetIdx()])
        pairs = list(itertools.combinations(symbols, 2))
        for pair in pairs:
           aa=f"{pair[0]}-{atom_assigned_types[atom_i.GetIdx()]}-{pair[1]}"
           par_aa=angle_params[aa]
           list_angle_params_to_print.append((aa,par_aa[0],par_aa[1])) 
    else:
        continue

print(f"Angles summary\n",list_angle_params_to_print)

#Improper dihedrals
list_impro_dihedrals_params_to_print=[]
for atom_i in mol.GetAtoms():
    neighbors=atom_i.GetNeighbors()
    if atom_i.GetNeighbors() == 3:
        symbols=[]
        for i in atom_i.GetNeighbors():
            symbols.append(atom_assigned_types[i.GetIdx()])
        tris = list(itertools.permutations(symbols, 3))
        for tri in tris:
            aa=f"{atom_assigned_types[atom_i.GetIdx()]}-{tri[0]}-{tri[1]}-{tri[2]}"
            if aa in improper_dihedral_params.keys():
                par_aa=improper_dihedral_params[aa]
                list_impro_dihedrals_params_to_print.append((aa,par_aa[0],par_aa[1],par_aa[2])) 
            else:
                continue
    else:
        continue
print(f"Improper dihedrals summary\n",list_impro_dihedrals_params_to_print)

#Proper dihedrals
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
                                    list_dihedrals_params_to_print.append((aa,par_aa[0],par_aa[1],par_aa[2]))
                                elif bb in dihedral_params.keys():
                                    par_bb=dihedral_params[bb]
                                    list_dihedrals_params_to_print.append((bb,par_bb[0],par_bb[1],par_bb[2]))
                                else:
                                    continue
                            else:
                                continue
                        else:
                            continue
                    else:
                        continue

print(f"Proper dihedrals summary\n",list_dihedrals_params_to_print)

#LJ params
list_LJ_params_to_print=[]
for atom_i in mol.GetAtoms():
    idx=atom_i.GetIdx()
    type_i=atom_assigned_types[idx]
    par_lj=LJ_params[type_i]
    list_LJ_params_to_print.append((idx,type_i,par_lj[0],par_lj[1]))

print(f"LJ params summary\n",list_LJ_params_to_print)
