from rdkit import Chem
import pandas as pd

def assign_atom_types(mol, smarts_to_type):
    """
    Assign atom types to a molecule based on SMARTS patterns.

    Parameters:
        mol: RDKit molecule
        smarts_to_type: dict mapping SMARTS -> atom type

    Returns:
        atom_types: list of atom types corresponding to atoms in mol
    """
    atom_types = [None] * mol.GetNumAtoms()

    for smarts, atom_type in smarts_to_type.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            raise ValueError(f"Invalid SMARTS: {smarts}")
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            for atom_idx in match:
                # Assign type only if not already assigned
                if mol.GetAtomWithIdx(atom_idx).GetSymbol()[0] == atom_type.capitalize()[0]: 
                    if atom_types[atom_idx] is None:
                        atom_types[atom_idx] = atom_type

    # Optionally, check for unassigned atoms
    for i, atype in enumerate(atom_types):
        if atype is None:
            atom_types[i] = "UNDEFINED"

    return atom_types

# Load molecule
mol = Chem.MolFromSmiles("CCO")  # ethanol
mol = Chem.AddHs(mol)
# Define SMARTS -> atom type
data=pd.read_csv("gaff.prm",sep='\s+',header=None,names=["atom","smarts","type","notes"])
data1=data["smarts"].to_list()
data2=data["type"].to_list()
smarts_to_type = { data1[i] : data2[i] for i in range(len(data2)) }
smarts_to_type_sorted = dict(
    sorted(smarts_to_type.items(), key=lambda x: -len(x[0]))
)
# Assign atom types
types = assign_atom_types(mol, smarts_to_type_sorted)

# Print
for atom, t in zip(mol.GetAtoms(), types):
    print(atom.GetSymbol(), t)

