# auto_typer.py (Version 10.1 - Fixed Guanidinium Detection)
import json
import argparse
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# --- FINALIZED RULES ---
TYPE_RULES = [
    # CHARGED SPECIES (Highest Precedence)
    ("N_PLUS", "[N+1]"),
    ("O_COO", "[O-1]"),

    # HYDROGENS
    ("HA", "[H](~c)"),
    ("H_N", "[H](~N)"),
    ("H_O", "[H](~O)"),
    ("H_C", "[H](~C)"),

    # SPECIFIC NEUTRAL HEAVY ATOMS
    ("CA", "[c]"),
    ("C_CO", "[C;D3](=[O;D1])"),
    ("N_AMINE", "[N;H2;D3]"),
    ("O_H", "[O;H2;D2]"),
    ("O_H", "[O;H1;D2]"),
    ("O_CO", "[O;D1]=[C]"),

    # GENERAL FALLBACKS
    ("CT", "[C;X4;!c]"),
    ("N_AMINE", "[N;X3;!N+]"),
]

def assign_special_groups(mol, assigned_types, smiles_string):
    """
    A dedicated function to handle complex functional groups that are poorly
    described by SMARTS. This runs before the main loop.
    """
    # --- Fixed Guanidinium Detection ---
    # Pattern for guanidinium: C with imine (C=N) + 2 amines, common in arginine
    guanidinium_pat = Chem.MolFromSmarts("[C](=[N])([N])[N]")
    if guanidinium_pat:
        matches = mol.GetSubstructMatches(guanidinium_pat)
        for match in matches:
            c_idx, n1_idx, n2_idx, n3_idx = match
            print(f"  [Fix] Found Guanidinium group: C={c_idx}, Ns=[{n1_idx},{n2_idx},{n3_idx}]")
            assigned_types[c_idx] = "C_GUA"
            assigned_types[n1_idx] = "N_GUA"
            assigned_types[n2_idx] = "N_GUA"
            assigned_types[n3_idx] = "N_GUA"

    # Fallback: any carbon with exactly 3 nitrogen neighbors
    if not guanidinium_pat or not mol.GetSubstructMatches(guanidinium_pat):
        for i, atom in enumerate(mol.GetAtoms()):
            if atom.GetSymbol() == 'C' and i not in assigned_types:
                n_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'N']
                if len(n_neighbors) == 3:
                    print(f"  [Fix] Found Guanidinium by neighbor count: C={i}, Ns={[n.GetIdx() for n in n_neighbors]}")
                    assigned_types[i] = "C_GUA"
                    for n in n_neighbors:
                        assigned_types[n.GetIdx()] = "N_GUA"

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
                c_atom = neighbor
                for c_neighbor in c_atom.GetNeighbors():
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

    # PASS 1: Handle complex, multi-atom functional groups first.
    assign_special_groups(mol, assigned_types, smiles_string)

    # PASS 2: Main rule-based assignment loop.
    for type_name, smarts in TYPE_RULES:
        pattern = Chem.MolFromSmarts(smarts)
        if not pattern: continue
        matches = mol.GetSubstructMatches(pattern)
        for match_tuple in matches:
            if not match_tuple: continue
            # Fixed: Handle single atom patterns correctly
            for idx in match_tuple:
                if idx not in assigned_types:
                    assigned_types[idx] = type_name
                    break  # Only assign to first unassigned atom in match

    # PASS 3: Post-processing fixes for resonance effects.
    fix_carboxylate_resonance(mol, assigned_types)

    # Final assignment to the JSON data.
    for i, atom_spec in enumerate(atoms_data):
        if i in assigned_types:
            atom_spec['type_name'] = assigned_types[i]
        else:
            # We add a print here to be notified of any future gaps in our rules.
            print(f"Warning: No rule matched atom {atom_spec['id']}. Check rules.", file=sys.stderr)
            atom_spec['type_name'] = atom_spec['element']

    bonds_data, double_bonds_data = [], []
    for bond in mol.GetBonds():
        a1_id = f"{mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()}{bond.GetBeginAtomIdx() + 1}"
        a2_id = f"{mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()}{bond.GetEndAtomIdx() + 1}"
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonds_data.append([a1_id, a2_id])
        else:
            bonds_data.append([a1_id, a2_id])

    return {"name": name, "atoms": atoms_data, "bonds": bonds_data, "double_bonds": double_bonds_data}

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
