# auto_typer.py (Version 12.0 - Final Amide/Ester/Acid Fix)
import json
import argparse
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# --- FINALIZED RULES ---
# This order is definitive.
TYPE_RULES = [
    # CHARGED SPECIES (Highest Precedence)
    ("N_PLUS", "[N+1]"),
    ("O_COO", "[O-1]"),

    # HYDROGENS
    ("HA", "[H](~c)"),
    ("H_N", "[H](~N)"),
    ("H_O", "[H](~O)"),
    ("H_C", "[H](~C)"),

    # SPECIAL HYBRIDIZATIONS
    ("C_NITRILE", "[C;D2]#[N;D1]"),
    ("N_NITRILE", "[N;D1]#[C;D2]"),

    #Carboxylic Acid Oxygen (must be before generic hydroxyls)
    ("O_ACID", "[O;H1;D2](C(=O))"),

    # SPECIFIC NEUTRAL HEAVY ATOMS - ORDER IS CRITICAL
    # Amide Nitrogen (must be before generic amines)
    ("N_AMIDE", "[N;D3](C(=O))"),
    # Ether/Ester Oxygen (must be before generic hydroxyls)
    ("O_ETHER", "[O;D2](C)C"),
    # Carbonyl Carbon
    ("C_CO", "[C;D3](=[O;D1])"),
    # Aromatic Carbon
    ("CA", "[c]"),

    # Generic Amines/Hydroxyls
    ("N_AMINE", "[N;H3;D3]"),     # Explicitly matches NH3
    ("N_AMINE", "[N;H2;D3]"),
    ("O_H", "[O;H2;D2]"),
    ("O_H", "[O;H1;D2]"),
    ("O_CO", "[O;D1]=[C]"),

    # GENERAL FALLBACKS (Lowest Precedence)
    ("CT", "[C;X4;!c]"),
    ("N_AMINE", "[N;X3;!N+]"),
]

def assign_special_groups(mol, assigned_types):
    """
    A dedicated function to handle complex functional groups procedurally.
    """
    # --- Guanidinium Group in Arginine ---
    # The procedural check for a carbon bonded to 3 nitrogens is the most robust.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIdx() not in assigned_types:
            neighbors = atom.GetNeighbors()
            n_neighbors = [n for n in neighbors if n.GetSymbol() == 'N']
            if len(n_neighbors) == 3:
                c_idx = atom.GetIdx()
                assigned_types[c_idx] = "C_GUA"
                for n_atom in n_neighbors:
                    assigned_types[n_atom.GetIdx()] = "N_GUA"
                break # Assume only one guanidinium group

    # --- Nitro Group ---
    nitro_pat = Chem.MolFromSmarts("[N+](=O)[O-]")
    if nitro_pat:
        matches = mol.GetSubstructMatches(nitro_pat)
        for match in matches:
            if len(match) == 3:
                n_idx, o1_idx, o2_idx = match
                assigned_types[n_idx], assigned_types[o1_idx], assigned_types[o2_idx] = "N_NO2", "O_NO2", "O_NO2"

def fix_carboxylate_resonance(mol, assigned_types):
    carboxylate_oxygens_indices = [i for i, t in assigned_types.items() if t == 'O_COO']
    for o_coo_idx in carboxylate_oxygens_indices:
        o_atom = mol.GetAtomWithIdx(o_coo_idx)
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                for c_neighbor in neighbor.GetNeighbors():
                    if c_neighbor.GetSymbol() == 'O' and c_neighbor.GetIdx() != o_coo_idx:
                        if assigned_types.get(c_neighbor.GetIdx()) == 'O_CO':
                            assigned_types[c_neighbor.GetIdx()] = 'O_COO'
                break

def build_molecule_from_smiles(smiles_string, name):
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol: return None

    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)

    atoms_data = []
    conformer = mol.GetConformer()
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        atom_id = f"{atom.GetSymbol()}{atom.GetIdx() + 1}"
        atoms_data.append({"id": atom_id, "type_name": "", "element": atom.GetSymbol(), "pos": [pos.x / 10.0, pos.y / 10.0, pos.z / 10.0]})

    assigned_types = {}

    assign_special_groups(mol, assigned_types)

    for type_name, smarts in TYPE_RULES:
        pattern = Chem.MolFromSmarts(smarts)
        if not pattern: continue
        matches = mol.GetSubstructMatches(pattern)
        for match_tuple in matches:
            # This loop now correctly handles single-atom patterns.
            idx = match_tuple[0]
            if idx not in assigned_types:
                assigned_types[idx] = type_name

    fix_carboxylate_resonance(mol, assigned_types)

    for i, atom_spec in enumerate(atoms_data):
        if i in assigned_types:
            atom_spec['type_name'] = assigned_types[i]
        else:
            print(f"Warning: No rule matched atom {atom_spec['id']}. Using element fallback.", file=sys.stderr)
            atom_spec['type_name'] = atom_spec['element']

    bonds_list = []
    for bond in mol.GetBonds():
        a1_id = f"{mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()}{bond.GetBeginAtomIdx() + 1}"
        a2_id = f"{mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()}{bond.GetEndAtomIdx() + 1}"
        bond_type = bond.GetBondType()
        order_str = "Single"
        if bond_type == Chem.BondType.DOUBLE: order_str = "Double"
        elif bond_type == Chem.BondType.TRIPLE: order_str = "Triple"
        bonds_list.append({"atoms": [a1_id, a2_id], "order": order_str})

    return {"name": name, "atoms": atoms_data, "bonds": bonds_list}

def main():
    parser = argparse.ArgumentParser(description="Generate a typed molecule JSON from a SMILES string.")
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--name", default="Generated Molecule")
    args = parser.parse_args()
    molecule_json = build_molecule_from_smiles(args.smiles, args.name)
    if molecule_json:
        print(json.dumps(molecule_json, indent=2))
    else:
        print(json.dumps({"error": f"Could not parse SMILES: {args.smiles}"}), file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
