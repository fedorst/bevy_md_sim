# auto_typer.py
import json
import argparse
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# THE FIX: Load the force field and create a lookup for defined bonds.
# This lookup will be our primary source of truth for bond orders.
FORCE_FIELD_PATH = os.path.join(os.path.dirname(__file__), 'assets', 'force_field.json')
with open(FORCE_FIELD_PATH) as f:
    FORCE_FIELD = json.load(f)

BOND_ORDER_MAP = {}
for bond_data in FORCE_FIELD['bonds']:
    types = tuple(sorted(bond_data['types']))
    order = bond_data.get('order', 'Single')
    BOND_ORDER_MAP[types] = order

# Use the original, proven TYPE_RULES
TYPE_RULES = [
    ("N_PLUS", "[N+1]"),("O_COO", "[O-1]"),("HA", "[H](~c)"),("H_N", "[H](~N)"),
    ("H_O", "[H](~O)"),("H_C", "[H](~C)"),("C_NITRILE", "[C;D2]#[N;D1]"),
    ("N_NITRILE", "[N;D1]#[C;D2]"),("O_ACID", "[O;H1;D2](C(=O))"),
    ("N_AMIDE", "[N;D3](C(=O))"),("O_ETHER", "[O;D2](C)C"),
    ("C_CO", "[C;D3](=[O;D1])"),("CA", "[c]"),("N_AMINE", "[N;H3;D3]"),
    ("N_AMINE", "[N;H2;D3]"),("O_H", "[O;H2;D2]"),("O_H", "[O;H1;D2]"),
    ("O_CO", "[O;D1]=[C]"),("CT", "[C;X4;!c]"),("N_AMINE", "[N;X3;!N+]"),
]

def manually_type_histidine_ring(mol, assigned_types):
    """
    Finds the 5-membered imidazole ring and manually types all of its atoms,
    including attached hydrogens. This is more robust than SMARTS for this case.
    """
    ring_info = mol.GetRingInfo()
    for ring_indices in ring_info.AtomRings():
        if len(ring_indices) == 5:
            # We found the 5-membered ring.
            for atom_idx in ring_indices:
                atom = mol.GetAtomWithIdx(atom_idx)
                # Manually type the heavy atoms in the ring.
                if atom.GetSymbol() == 'N':
                    assigned_types[atom_idx] = 'N_IM'
                elif atom.GetSymbol() == 'C':
                    assigned_types[atom_idx] = 'C_IM'

                # CRUCIAL: Also manually type their attached hydrogens now.
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'H':
                        # The general H_N/H_C rules are fine for the label.
                        if atom.GetSymbol() == 'N':
                            assigned_types[neighbor.GetIdx()] = 'H_N'
                        elif atom.GetSymbol() == 'C':
                            assigned_types[neighbor.GetIdx()] = 'H_C'
            return # Stop after finding the first 5-membered ring.

def assign_special_groups(mol, assigned_types):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIdx() not in assigned_types:
            neighbors = atom.GetNeighbors()
            if len([n for n in neighbors if n.GetSymbol() == 'N']) == 3:
                assigned_types[atom.GetIdx()] = "C_GUA"
                for n_atom in neighbors:
                    if n_atom.GetSymbol() == 'N': assigned_types[n_atom.GetIdx()] = "N_GUA"
                break
    nitro_pat = Chem.MolFromSmarts("[N+](=O)[O-]")
    if nitro_pat:
        for match in mol.GetSubstructMatches(nitro_pat):
            if len(match) == 3:
                n_idx, o1_idx, o2_idx = match
                assigned_types[n_idx], assigned_types[o1_idx], assigned_types[o2_idx] = "N_NO2", "O_NO2", "O_NO2"

def fix_carboxylate_resonance(mol, assigned_types):
    for o_coo_idx in [i for i, t in assigned_types.items() if t == 'O_COO']:
        for neighbor in mol.GetAtomWithIdx(o_coo_idx).GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                for c_neighbor in neighbor.GetNeighbors():
                    if c_neighbor.GetSymbol() == 'O' and c_neighbor.GetIdx() != o_coo_idx and assigned_types.get(c_neighbor.GetIdx()) == 'O_CO':
                        assigned_types[c_neighbor.GetIdx()] = 'O_COO'
                break

def build_molecule_from_smiles(smiles_string, name):
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol: return None

    mol = Chem.AddHs(mol)

    # THE FIX: Perform atom typing BEFORE sanitization and embedding.
    atoms_data = [{"id": f"{a.GetSymbol()}{a.GetIdx() + 1}", "type_name": "", "element": a.GetSymbol(), "pos": [0,0,0]} for a in mol.GetAtoms()]
    assigned_types = {}
    assign_special_groups(mol, assigned_types)

    if "C1=CNC=N1" in smiles_string:
        manually_type_histidine_ring(mol, assigned_types)


    for type_name, smarts in TYPE_RULES:
        pattern = Chem.MolFromSmarts(smarts)
        if not pattern: continue
        for match_tuple in mol.GetSubstructMatches(pattern):
            if match_tuple[0] not in assigned_types:
                assigned_types[match_tuple[0]] = type_name


    fix_carboxylate_resonance(mol, assigned_types)

    atom_type_map_by_index = {}
    for i, atom_spec in enumerate(atoms_data):
        atom_spec['type_name'] = assigned_types.get(i, atom_spec['element'])
        atom_type_map_by_index[i] = atom_spec['type_name']

    # Now that we have types, we can assign bond orders from our force field.
    bonds_list = []
    for bond in mol.GetBonds():
        a1_idx, a2_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        type1, type2 = atom_type_map_by_index[a1_idx], atom_type_map_by_index[a2_idx]
        bond_key = tuple(sorted((type1, type2)))
        order_str = BOND_ORDER_MAP.get(bond_key, "Single")

        a1_id = f"{mol.GetAtomWithIdx(a1_idx).GetSymbol()}{a1_idx + 1}"
        a2_id = f"{mol.GetAtomWithIdx(a2_idx).GetSymbol()}{a2_idx + 1}"
        bonds_list.append({"atoms": [a1_id, a2_id], "order": order_str})

    # Now that typing is done, we can sanitize and generate coordinates.
    try: Chem.SanitizeMol(mol)
    except Exception: pass
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)

    conformer = mol.GetConformer()
    for i, atom_spec in enumerate(atoms_data):
        pos = conformer.GetAtomPosition(i)
        atom_spec['pos'] = [pos.x / 10.0, pos.y / 10.0, pos.z / 10.0]

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
