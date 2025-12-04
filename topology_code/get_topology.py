"""Atom Typer and Topology File Generator
Starting from SMILES (python get_topology.py -s/--smi'SMILES' -o/--output output_filename) or
PDB file (python get_topology.py -p/--pdb coord_file.pdb -o/--output output_filename), it assigns the atom types
according to predefined SMARTS patterns. Then, it proceeds to scan all the molecule bonds,
angles, proper and improper dihedral to assign the corresponding parameters to each and store them in lists. 
It also assigns parameteres for Lennard-Jones potential and it calculates the RESP charges for each
atom at HF/STO-3G level using psi4 and its resp plugin. Constrains in the RESP calculation may be given with the input
flag -c/--constraints Finally, the whole set of parameters needed for the 
chosen molecule are printed in a *.top file, ordered by type of interaction, and the XYZ coordinates are saved
into *.xyz file. 

The atom types are taken from GAFF.

In the output file, atoms' index count starts at 1, instead of 0, to improve compatibility with Fortran scripts.
Example: python get_topology.py --smi 'c1ccccc1' -o benzene

Author: Leonardo Di Ciano (2025)"""

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from force_field_data import mass_dict, bond_params, angle_params, improper_dihedral_params, dihedral_params, LJ_params
import sys
import itertools
import warnings
import psi4
import resp
import os
import glob
import argparse
import ast

parser=argparse.ArgumentParser(prog="Topology file generator")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-s","--smi", type=str, help="Input SMILES string")
group.add_argument("-p","--pdb", type=str, help="Input PDB file path")
parser.add_argument('-o','--output',type=str,metavar='filename',help=f"Give the name for topology and XYZ coordinate file")
parser.add_argument('-c','--constraints',type=str,metavar=None,help=f"Give a list of constrained groups of atoms for RESP charges calculation, i.e. '[[1,2,][3,4,5,6,7,8]]' for ethane")
args=parser.parse_args()

constraints = None
if args.constraints:
    try:
        constraints = ast.literal_eval(args.constraints)
    except:
        pass
     
#Avoid pandas FutureWarnings
warnings.filterwarnings('ignore',category=FutureWarning)
warnings.filterwarnings('ignore',category=SyntaxWarning)

#SMILES or XYZ file run mode
if args.smi :
    mol=Chem.MolFromSmiles(sys.argv[2])
    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
elif args.pdb:
    mol=Chem.MolFromPDBFile(sys.argv[2])
    AllChem.SanitizeMol(mol)
    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
else:
    raise KeyError

"""Atom Types Definition
The assignment is generally based on combinations of type of neighbor atoms, 
atoms hybridization and aromatic nature."""
atom_assigned_types={} #dict with index and atom types
list_atom_types_to_print = [None] * mol.GetNumAtoms()

typer_data=pd.read_csv("FF_csv/atom_types.ff",sep='\s+',header=None,names=["atom","smarts","type","notes"])
data1=typer_data["smarts"].to_list()
data2=typer_data["type"].to_list()
smarts_to_type = { data1[i] : data2[i] for i in range(len(data2)) }
smarts_to_type_sorted = dict(
    sorted(smarts_to_type.items(), key=lambda x: -len(x[0]))
)
for smarts, atom_type in smarts_to_type_sorted.items():
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f"Invalid SMARTS: {smarts}")
    matches = mol.GetSubstructMatches(pattern)
    for match in matches:
        for atom_idx in match:
            symb=mol.GetAtomWithIdx(atom_idx).GetSymbol()
            if symb[0] == atom_type.capitalize()[0]: #Check that is the correct element 
                if list_atom_types_to_print[atom_idx] is None:
                    list_atom_types_to_print[atom_idx] = [atom_idx + 1,symb, atom_type, float(mass_dict[symb])]
                    atom_assigned_types[atom_idx]=atom_type



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
            symbols.append(i.GetIdx())
        pairs = list(itertools.combinations(symbols, 2))
        for pair in pairs:
           aa=f"{atom_assigned_types[pair[0]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[pair[1]]}"
           par_aa=angle_params[aa]
           list_angle_params_to_print.append((aa,pair[0]+1,atom_i.GetIdx()+1,pair[1]+1,par_aa[0],par_aa[1])) 
    else:
        continue

print(f"Angles summary\n",list_angle_params_to_print)

"""Improper Diehedrals assignments
Assigned by iterating over the possible combinations of a central atom and three *different* neighbors,
itertools.permutations allows repetitions, in order to get a correctly ordered set to grep from the 
force_field_data library"""
quartet_list=[]
list_impro_dihedrals_params_to_print=[]
for atom_i in mol.GetAtoms():
    neighbors=atom_i.GetNeighbors()
    if len(atom_i.GetNeighbors()) == 3:
        symbols=[]
        for i in atom_i.GetNeighbors():
            symbols.append(i.GetIdx())
        tris = list(itertools.permutations(symbols, 3))
        for tri in tris:
            aa=f"{atom_assigned_types[tri[0]]}-{atom_assigned_types[tri[1]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[2]]}"
            bb=f"X-{atom_assigned_types[tri[0]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[1]]}"
            cc=f"X-{atom_assigned_types[tri[1]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[2]]}"
            dd=f"X-{atom_assigned_types[tri[0]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[2]]}"
            ee=f"X-X-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[0]]}"
            ff=f"X-X-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[1]]}"
            gg=f"X-X-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[2]]}"
            if aa in improper_dihedral_params.keys():
                par_aa=improper_dihedral_params[aa]
                list_impro_dihedrals_params_to_print.append((aa,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_aa[0],par_aa[1],par_aa[2],par_aa[3])) 
                quartet_list.extend(itertools.permutations([tri[0], tri[1], atom_i.GetIdx(), tri[2]], 4))
            elif bb in improper_dihedral_params.keys():
                par_bb=improper_dihedral_params[bb]
                list_impro_dihedrals_params_to_print.append((bb,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_bb[0],par_bb[1],par_bb[2],par_bb[3])) 
            elif cc in improper_dihedral_params.keys():
                par_cc=improper_dihedral_params[cc]
                list_impro_dihedrals_params_to_print.append((cc,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_cc[0],par_cc[1],par_cc[2],par_cc[3])) 
            elif dd in improper_dihedral_params.keys():
                par_dd=improper_dihedral_params[dd]
                list_impro_dihedrals_params_to_print.append((dd,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_dd[0],par_dd[1],par_dd[2],par_dd[3])) 
            elif ee in improper_dihedral_params.keys():
                par_ee=improper_dihedral_params[ee]
                list_impro_dihedrals_params_to_print.append((ee,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_ee[0],par_ee[1],par_ee[2],par_ee[3])) 
            elif ff in improper_dihedral_params.keys():
                par_ff=improper_dihedral_params[ff]
                list_impro_dihedrals_params_to_print.append((ff,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_ff[0],par_ff[1],par_ff[2],par_ff[3])) 
            elif gg in improper_dihedral_params.keys():
                par_gg=improper_dihedral_params[gg]
                list_impro_dihedrals_params_to_print.append((gg,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_gg[0],par_gg[1],par_gg[2],par_gg[3])) 
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
quartet_list=[]
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
                                cc=f"X-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[atom_j.GetIdx()]}-X"
                                if not [i_nghb.GetIdx(),atom_i.GetIdx(),atom_j.GetIdx(),j_nghb.GetIdx()] in quartet_list:
                                    if aa in dihedral_params.keys() :
                                        params_aa=dihedral_params[aa]
                                        quartet_list.append([i_nghb.GetIdx(),atom_i.GetIdx(),atom_j.GetIdx(),j_nghb.GetIdx()])
                                        quartet_list.append([j_nghb.GetIdx(),atom_j.GetIdx(),atom_i.GetIdx(),i_nghb.GetIdx()])
                                        for par_aa in params_aa:
                                            list_dihedrals_params_to_print.append((aa,i_nghb.GetIdx()+1,atom_i.GetIdx()+1,atom_j.GetIdx()+1,j_nghb.GetIdx()+1,par_aa[0],par_aa[1],par_aa[2],par_aa[3]))
                                    elif bb in dihedral_params.keys() :
                                        params_bb=dihedral_params[bb]
                                        quartet_list.append([i_nghb.GetIdx(),atom_i.GetIdx(),atom_j.GetIdx(),j_nghb.GetIdx()])
                                        quartet_list.append([j_nghb.GetIdx(),atom_j.GetIdx(),atom_i.GetIdx(),i_nghb.GetIdx()])
                                        for par_bb in params_bb:
                                            list_dihedrals_params_to_print.append((bb,j_nghb.GetIdx()+1,atom_j.GetIdx()+1,atom_i.GetIdx()+1,i_nghb.GetIdx()+1,par_bb[0],par_bb[1],par_bb[2],par_bb[3]))
                                    elif cc in dihedral_params.keys():
                                        params_cc=dihedral_params[cc]
                                        quartet_list.append([i_nghb.GetIdx(),atom_i.GetIdx(),atom_j.GetIdx(),j_nghb.GetIdx()])
                                        quartet_list.append([j_nghb.GetIdx(),atom_j.GetIdx(),atom_i.GetIdx(),i_nghb.GetIdx()])
                                        for par_cc in params_cc:
                                            list_dihedrals_params_to_print.append((cc,j_nghb.GetIdx()+1,atom_j.GetIdx()+1,atom_i.GetIdx()+1,i_nghb.GetIdx()+1,par_cc[0],par_cc[1],par_cc[2],par_cc[3]))
                                    else:
                                        continue
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
psi4.core.set_output_file('output.dat', False)
options = {'VDW_SCALE_FACTORS' : [1.4, 1.6, 1.8, 2.0],
           'VDW_POINT_DENSITY'  : 1.0,
           'RESP_A'             : 0.0005,
           'RESP_B'             : 0.1,
           'RESTRAINT'          : True,
           'IHFREE'             : False,
           'WEIGHT'             : [1, 1, 1],
           'BASIS_ESP'          : "STO-3G"
           }
# Call for first stage fit
charges1 = resp.resp([psi_mol], options)
# Change the value of the RESP parameter A
options['RESP_A'] = 0.0001
# Add constraint for atoms fixed in second stage fit
constraint_charge = []
for i in range(len(charges1)):
    constraint_charge.append([charges1[1][i], [i+1]])
options['constraint_charge'] = constraint_charge
if constraints:
    options['constraint_group'] = constraints
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
if args.output:
    file=f"{args.output}.top"
else:
    file="topology.top"
with open(file,"w+") as f:
    f.write(f"[AtomTypes] {len(list_atom_types_to_print)}\n")
    types_df=pd.DataFrame(list_atom_types_to_print,columns=["Idx","Symb","Type","MW"],index=None)
    types_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[Bonds] {len(list_bond_params_to_print)}\n")
    bonds_df=pd.DataFrame(list_bond_params_to_print,columns=["Idx1","Idx2","b0","kb"],index=None)
    bonds_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[Angles] {len(list_angle_params_to_print)}\n")
    angles_df=pd.DataFrame(list_angle_params_to_print,columns=["Types","Idx1","Idx2","Idx3","th0","cth"],index=None)
    angles_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[ImproperDihedrals] {len(list_impro_dihedrals_params_to_print)}\n")
    impdie_df=pd.DataFrame(list_impro_dihedrals_params_to_print,columns=["Types","Idx1","Idx2","Idx3","Idx4","div","phase","kd","pn"],index=None)
    impdie_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[Dihedrals] {len(list_dihedrals_params_to_print)}\n")
    die_df=pd.DataFrame(list_dihedrals_params_to_print,columns=["Types","Idx1","Idx2","Idx3","Idx4","div","phase","kd","pn"],index=None)
    die_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[LJ]\n")
    lj_df=pd.DataFrame(list_LJ_params_to_print,columns=["Idx","sigma","epsilon"],index=None)
    lj_df.to_csv(f,sep="\t",header=False,index=False,float_format='%.5f',escapechar=None)
    f.write(f"\n[Charges]\n")
    chg_df=pd.DataFrame(list_resp_charges_to_print,columns=["Idx","charge"],index=None)
    chg_df.to_csv(f,sep="\t",header=False,index=False,float_format='%.3f',escapechar=None)
    f.write(f" ")
print(f"Topology file written in {file}")



if args.output:
    filename=f"{args.output}.xyz"
else:
    filename="coord.xyz"
Chem.MolToXYZFile(mol,filename)
print(f"XYZ coordinates written in {filename}")


