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

The atom types are taken from GAFF Version 1.4, March 2010.

In the output file, atoms' index count starts at 1, instead of 0, to improve compatibility with Fortran scripts.
Example: python get_topology.py --smi 'c1ccccc1' -o benzene

Author: Leonardo Di Ciano (2025)"""

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from force_field_data import mass_dict, bond_params, angle_params, improper_dihedral_params, dihedral_params, LJ_params, smarts_to_type
import itertools
import warnings
import psi4
import resp
import os
import glob
import argparse
import ast
from get_mulliken_charges import Mulliken_Charges

parser=argparse.ArgumentParser(prog="get_topology.py")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-s","--smi", type=str, help="Input SMILES string")
group.add_argument("-p","--pdb", type=str, help="Input PDB file path")
parser.add_argument('-o','--output',type=str,metavar='filename',help=f"Give the name for topology and XYZ coordinate file")
group2 = parser.add_mutually_exclusive_group(required=False)
group2.add_argument('-m','--mulliken',action='store_true',help=f"Select Mulliken Charges calculation to obtain partial charges")
group2.add_argument('-r','--resp',action='store_true',help=f"Select RESP charges calculation to obtain partial charges")
group2.add_argument('-a','--assign',type=str,metavar=None,help=f"Give a list of atom types and corresponding charge to assign")
parser.add_argument('-c','--constraints',type=str,metavar=None,help=f"Give a list of constrained groups of atoms for RESP charges calculation, i.e. '[[1,2,][3,4,5,6,7,8]]' for ethane")
parser.add_argument('-d','--debug',action='store_true',help=f"Activate debug level of printing")
args=parser.parse_args()

if not args.mulliken and not args.assign:
    args.resp = True
if args.constraints and not args.resp:
    parser.error("-c/--constraints can only be used together with -r/--resp")
elif args.constraints and args.resp:
    try:
        constraints = ast.literal_eval(args.constraints)
        print(constraints)
    except:
        constraints = None
else:
    constraints = None

if args.assign:
    try:
        assignments = ast.literal_eval(args.assign)
        print(assignments)
    except:
        assignments = None
else:
    assignments = None
     
#Avoid pandas FutureWarnings
warnings.filterwarnings('ignore',category=FutureWarning)
warnings.filterwarnings('ignore',category=SyntaxWarning)

#SMILES or XYZ file run mode
if args.smi :
    mol=Chem.MolFromSmiles(args.smi)
    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
elif args.pdb:
    mol=Chem.MolFromPDBFile(args.pdb,removeHs=False)
    AllChem.SanitizeMol(mol)
else:
    raise KeyError


if args.output:
    xyzfilename=f"{args.output}.xyz"
else:
    xyzfilename="coord.xyz"
Chem.MolToXYZFile(mol,xyzfilename)
print(f"XYZ coordinates written in {xyzfilename}")


"""Atom Types Definition
The assignment is generally based on combinations of type of neighbor atoms, 
atoms hybridization and aromatic nature."""
atom_assigned_types={} #dict with index and atom types
list_atom_types_to_print = [None] * mol.GetNumAtoms()

for smarts, atom_type in smarts_to_type.items():
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

for i in range(len(list_atom_types_to_print)):     
    if list_atom_types_to_print[i] == None:              
        if mol.GetAtomWithIdx(i).GetSymbol() == "H":
            symbols=[]
            for j in mol.GetAtomWithIdx(i).GetNeighbors():
                symbols.append(j.GetSymbol())
            if "O" in symbols:
                list_atom_types_to_print[i] = [i + 1,"H", "ho", float(mass_dict["H"])]
                atom_assigned_types[i] = "ho"   
            elif "N" in symbols:
                list_atom_types_to_print[i] = [i + 1,"H", "hn", float(mass_dict["H"])]
                atom_assigned_types[i] = "hn"
            else:
                list_atom_types_to_print[i] = [i + 1,"H", "hc", float(mass_dict["H"])]
                atom_assigned_types[i] = "hc"   

print(f"Atom Types assigned\n")
print_df=pd.DataFrame(list_atom_types_to_print,columns=["Index","Element","Type","Mass (amu)"])
print(print_df.to_csv(sep="\t",index=False))

"""Bonds assignments
Assigned simply by iterating over the bonds from rdkit and getting corresponding params
from the force_field_data library"""
list_bond_params_to_print=[]
for bond_i in mol.GetBonds():
    a1=bond_i.GetBeginAtomIdx()
    a2=bond_i.GetEndAtomIdx()
    if f"{atom_assigned_types[a1]}-{atom_assigned_types[a2]}" in bond_params.keys():
        par_12=bond_params[f"{atom_assigned_types[a1]}-{atom_assigned_types[a2]}"]
        list_bond_params_to_print.append((a1+1,a2+1,par_12[0],par_12[1]))
    elif f"{atom_assigned_types[a2]}-{atom_assigned_types[a1]}" in bond_params.keys():
        par_21=bond_params[f"{atom_assigned_types[a2]}-{atom_assigned_types[a1]}"]
        list_bond_params_to_print.append((a2+1,a1+1,par_21[0],par_21[1]))
    else:
        raise Exception(f"{atom_assigned_types[a1]}-{atom_assigned_types[a2]} bond not found")

if args.debug:
    print(f"Bond parameters assigned\n")
    print("k_b is the bond force constant in kcal/mol/Å^2")
    print(f"d_eq is the equilibrium bond distance in Å\n")
    print_df=pd.DataFrame(list_bond_params_to_print,columns=["Atom 1","Atom 2","k_b","d_eq"])
    print(print_df.to_csv(sep="\t",index=False,))
else:
    print(f"Bond parameters assigned: {len(list_bond_params_to_print)}")

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
            if aa in angle_params.keys():
                par_aa=angle_params[aa]
                list_angle_params_to_print.append((aa,pair[0]+1,atom_i.GetIdx()+1,pair[1]+1,par_aa[0],par_aa[1]))
            else:
                continue
    else:
        continue

if args.debug:
    print(f"Angle parameters assigned\n")
    print("k_theta is the angle force constant in kcal/mol/rad^2")
    print(f"theta_eq is the equilibrium angle in deg\n")
    print_df=pd.DataFrame(list_angle_params_to_print,columns=["Type","Atom 1","Atom 2","Atom 3","k_theta","theta_eq"])
    print(print_df.to_csv(sep="\t",index=False,))
else:
    print(f"Angle parameters assigned: {len(list_angle_params_to_print)}")


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


if args.debug:
    print(f"Proper dihedrals parameters assigned\n")
    print("k_phi is the dihedral energy barrier in kcal/mol")
    print("div is the divider for the energy barrier")
    print("phase is the phase angle in deg")
    print(f"n is the periodicity\n")
    print_df=pd.DataFrame(list_dihedrals_params_to_print,columns=["Type","Atom 1","Atom 2","Atom 3","Atom 4","divider","k_phi","phase","n"])
    print(print_df.to_csv(sep="\t",index=False,))
else:
    print(f"Proper dihedrals parameters assigned: {len(list_dihedrals_params_to_print)}")

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
            if (tri[0],tri[1],atom_i.GetIdx(),tri[2]) not in quartet_list:
                aa=f"{atom_assigned_types[tri[0]]}-{atom_assigned_types[tri[1]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[2]]}"
                bb=f"X-{atom_assigned_types[tri[0]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[1]]}"
                cc=f"X-{atom_assigned_types[tri[1]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[2]]}"
                dd=f"X-{atom_assigned_types[tri[0]]}-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[2]]}"
                ee=f"X-X-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[0]]}"
                ff=f"X-X-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[1]]}"
                gg=f"X-X-{atom_assigned_types[atom_i.GetIdx()]}-{atom_assigned_types[tri[2]]}"
                if aa in improper_dihedral_params.keys():
                    par_aa=improper_dihedral_params[aa]
                    list_impro_dihedrals_params_to_print.append((aa,tri[0]+1,tri[1]+1,atom_i.GetIdx()+1,tri[2]+1,par_aa[0],par_aa[1],par_aa[2])) 
                    quartet_list.extend(itertools.permutations([tri[0], tri[1], atom_i.GetIdx(), tri[2]], 4))
                elif bb in improper_dihedral_params.keys():
                    par_bb=improper_dihedral_params[bb]
                    list_impro_dihedrals_params_to_print.append((bb,tri[2]+1,tri[0]+1,atom_i.GetIdx()+1,tri[1]+1,par_bb[0],par_bb[1],par_bb[2])) 
                    quartet_list.extend(itertools.permutations([tri[0], tri[1], atom_i.GetIdx(), tri[2]], 4))
                elif cc in improper_dihedral_params.keys():
                    par_cc=improper_dihedral_params[cc]
                    list_impro_dihedrals_params_to_print.append((cc,tri[0]+1,tri[1]+1,atom_i.GetIdx()+1,tri[2]+1,par_cc[0],par_cc[1],par_cc[2])) 
                    quartet_list.extend(itertools.permutations([tri[0], tri[1], atom_i.GetIdx(), tri[2]], 4))
                elif dd in improper_dihedral_params.keys():
                    par_dd=improper_dihedral_params[dd]
                    list_impro_dihedrals_params_to_print.append((dd,tri[1]+1,tri[0]+1,atom_i.GetIdx()+1,tri[2]+1,par_dd[0],par_dd[1],par_dd[2])) 
                    quartet_list.extend(itertools.permutations([tri[0], tri[1], atom_i.GetIdx(), tri[2]], 4))
                elif ee in improper_dihedral_params.keys():
                    par_ee=improper_dihedral_params[ee]
                    list_impro_dihedrals_params_to_print.append((ee,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_ee[0],par_ee[1],par_ee[2])) 
                    quartet_list.extend(itertools.permutations([tri[0], tri[1], atom_i.GetIdx(), tri[2]], 4))
                elif ff in improper_dihedral_params.keys():
                    par_ff=improper_dihedral_params[ff]
                    list_impro_dihedrals_params_to_print.append((ff,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_ff[0],par_ff[1],par_ff[2]))
                    quartet_list.extend(itertools.permutations([tri[0], tri[1], atom_i.GetIdx(), tri[2]], 4))
                elif gg in improper_dihedral_params.keys():
                    par_gg=improper_dihedral_params[gg]
                    list_impro_dihedrals_params_to_print.append((gg,atom_i.GetIdx()+1,tri[0]+1,tri[1]+1,tri[2]+1,par_gg[0],par_gg[1],par_gg[2])) 
                    quartet_list.extend(itertools.permutations([tri[0], tri[1], atom_i.GetIdx(), tri[2]], 4))
                else:
                    continue
            else:
                continue
    else:
        continue

if args.debug:
    print(f"Improper dihedrals parameters assigned\n")
    print("k_phi is the dihedral energy barrier in kcal/mol")
    print("div is the divider for the energy barrier")
    print("phase is the phase angle in deg")
    print(f"n is the periodicity\n")
    print_df=pd.DataFrame(list_impro_dihedrals_params_to_print,columns=["Type","Atom 1","Atom 2","Atom 3","Atom 4","divider","k_phi","phase","n"])
    print(print_df.to_csv(sep="\t",index=False,))
else:
    print(f"Improper dihedrals parameters assigned: {len(list_impro_dihedrals_params_to_print)}")



"""Lennard-Jones assignments
The parameters are simply assigned for each atom type."""
list_LJ_params_to_print=[]
for atom_i in mol.GetAtoms():
    idx=atom_i.GetIdx()
    type_i=atom_assigned_types[idx]
    par_lj=LJ_params[type_i]
    list_LJ_params_to_print.append((idx+1,par_lj[0],par_lj[1]))

if args.debug:
    print(f"LJ parameters assigned\n")
    print("epsilon is the depth of the energy well in kcal/mol")
    print(f"sigma is the equilibrium distance in Å\n")
    print_df=pd.DataFrame(list_LJ_params_to_print,columns=["Index","epsilon","sigma"])
    print(print_df.to_csv(sep="\t",index=False,))
else:
    print(f"LJ parameters assigned: {len(list_LJ_params_to_print)}")

list_charges_to_print=[]
if args.resp:
    """RESP charges assignments
    The charges for Coulomb electrostatic interaction are not available in the library, but are calculated as 
    RESP charges. In this case, it uses the resp plugin of psi4 and there charges are obtained with two iterations,
    following https://github.com/cdsgroup/resp/blob/master/examples/example1.py """
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
            'BASIS_ESP'          : "6-31G*"
            }
    # Call for first stage fit
    charges1 = resp.resp([psi_mol], options)
    # Change the value of the RESP parameter A
    options['RESP_A'] = 0.0001
    # Add constraint for atoms fixed in second stage fit
    constraint_charge = []
    for i in range(len(charges1[1])):
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
        list_charges_to_print.append((i+1,charges2[1][i]))

    for i in glob.glob("*.dat"):
        os.remove(i)
    os.remove("results.out")

    if args.debug:
        print(f"RESP charges assigned\n")

        print_df=pd.DataFrame(list_charges_to_print,columns=["Index","Charge"])
        print(print_df.to_csv(sep="\t",index=False,))
    else:
        print(f"RESP charges assigned: {len(list_charges_to_print)}")
elif args.mulliken:
    """Mulliken charges assignment"""
    print("Calculating Mulliken Charges with Psi4")

    input=Mulliken_Charges(Chem.MolToXYZBlock(mol))
    charges=input.main()
    for i in range(len(charges)):
        list_charges_to_print.append((i+1,charges[i]))

    os.remove("PSI4_output.dat")

    if args.debug:
        print(f"Mulliken charges assigned\n")

        print_df=pd.DataFrame(list_charges_to_print,columns=["Index","Charge"])
        print(print_df.to_csv(sep="\t",index=False,))
    else:
        print(f"Mulliken charges assigned: {len(list_charges_to_print)}")
elif args.assign and assignments:
    """Charge assignment from input"""
    for i in range(len(list_atom_types_to_print)):
        list_charges_to_print.append((i+1,assignments[i%len(assignments)]))
    
    if args.debug:
        print(f"Charges assigned\n")

        print_df=pd.DataFrame(list_charges_to_print,columns=["Index","Charge"])
        print(print_df.to_csv(sep="\t",index=False,))
    else:
        print(f"Charges assigned: {len(list_charges_to_print)}")

"""The parameters are assembled with all the stored components and
written in a *.top file"""
if args.output:
    file=f"{args.output}.top"
else:
    file="topology.top"
with open(file,"w+") as f:
    f.write(f"[AtomTypes] {len(list_atom_types_to_print)}\n")
    types_df=pd.DataFrame(list_atom_types_to_print,index=None)
    types_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[Bonds] {len(list_bond_params_to_print)}\n")
    bonds_df=pd.DataFrame(list_bond_params_to_print,index=None)
    bonds_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[Angles] {len(list_angle_params_to_print)}\n")
    angles_df=pd.DataFrame(list_angle_params_to_print,index=None)
    angles_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[Dihedrals] {len(list_dihedrals_params_to_print)}\n")
    die_df=pd.DataFrame(list_dihedrals_params_to_print,index=None)
    die_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[ImproperDihedrals] {len(list_impro_dihedrals_params_to_print)}\n")
    impdie_df=pd.DataFrame(list_impro_dihedrals_params_to_print,index=None)
    impdie_df.to_csv(f,sep="\t",header=False,index=False,escapechar=None)
    f.write(f"\n[LJ]\n")
    lj_df=pd.DataFrame(list_LJ_params_to_print,index=None)
    lj_df.to_csv(f,sep="\t",header=False,index=False,float_format='%.5f',escapechar=None)
    f.write(f"\n[Charges]\n")
    chg_df=pd.DataFrame(list_charges_to_print,index=None)
    chg_df.to_csv(f,sep="\t",header=False,index=False,float_format='%.3f',escapechar=None)
    f.write(f" ")
print(f"Topology file written in {file}")



