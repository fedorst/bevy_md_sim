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
    ("N_PLUS", "[N+1]"),
    ("O_COO", "[O-1]"),
    ("N_AR", "[n]"),
    ("CA", "[c]"),
    ("S_THIOL", "[S;H1,H0;X2](-C)"),
    ("H_S", "[H](~S)"),
    ("S_THIOETHER", "[S;D2](-C)(-C)"),
    ("HA", "[H](~a)"),
    ("H_N", "[H](~N)"),
    ("H_O", "[H](~O)"),
    ("H_C", "[H](~C)"),
    ("C_NITRILE", "[C;D2]#[N;D1]"),
    ("N_NITRILE", "[N;D1]#[C;D2]"),
    ("O_ACID", "[O;H1;D2](C(=O))"),
    ("N_AMIDE", "[N;D3](C(=O))"),
    ("O_ETHER", "[O;D2](C)C"),
    ("C_CO", "[C;D3](=[O;D1])"),
    ("N_AMINE", "[N;H3;D3]"),
    ("N_AMINE", "[N;H2;D3]"),
    ("O_H", "[O;H2;D2]"),
    ("O_H", "[O;H1;D2]"),
    ("O_CO", "[O;D1]=[C]"),
    ("CT", "[C;X4;!c]"),
    ("N_AMINE", "[N;X3;!N+]"),
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


def rename_backbone_atom_ids(mol, atoms_data):
    """
    Finds the N-CA-C backbone and assigns standard IDs ("N", "CA", "C")
    to the 'id' field for easy lookup.
    """
    # A robust pattern for the amino acid backbone
    backbone_pattern = Chem.MolFromSmarts("[N;!H0;v3,v4]-[C;X4](-[H])-[C;X3](=[O;D1])")
    matches = mol.GetSubstructMatches(backbone_pattern)

    if not matches: # Fallback for Proline
        backbone_pattern = Chem.MolFromSmarts("[N;H1,H2;r5]-[C;X4;r5](-[C;X3](=[O;D1]))")
        matches = mol.GetSubstructMatches(backbone_pattern)
        if not matches: return

    n_idx, ca_idx, _, c_idx = matches[0][:4]

    for atom_spec in atoms_data:
        idx = atom_spec['idx']
        if idx == n_idx: atom_spec['id'] = 'N'
        elif idx == ca_idx: atom_spec['id'] = 'CA'
        elif idx == c_idx: atom_spec['id'] = 'C'

def build_molecule_from_smiles(smiles_string, name):
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol: return None

    mol = Chem.AddHs(mol)
    try: Chem.SanitizeMol(mol)
    except Exception: pass

    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)

    atoms_data = []
    conformer = mol.GetConformer()
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        atoms_data.append({
            "idx": atom.GetIdx(),
            "id": f"{atom.GetSymbol()}{atom.GetIdx() + 1}",
            "type_name": "", "element": atom.GetSymbol(),
            "pos": [pos.x / 10.0, pos.y / 10.0, pos.z / 10.0]
        })

    assigned_types = {}
    assign_special_groups(mol, assigned_types)
    if "C1=CNC=N1" in smiles_string: manually_type_histidine_ring(mol, assigned_types)
    for type_name, smarts in TYPE_RULES:
        pattern = Chem.MolFromSmarts(smarts)
        if not pattern: continue
        for match_tuple in mol.GetSubstructMatches(pattern):
            if match_tuple[0] not in assigned_types:
                assigned_types[match_tuple[0]] = type_name
    fix_carboxylate_resonance(mol, assigned_types)

    atom_type_map_by_index = {spec['idx']: assigned_types.get(spec['idx'], spec['element']) for spec in atoms_data}
    for spec in atoms_data:
        spec['type_name'] = atom_type_map_by_index[spec['idx']]

    rename_backbone_atom_ids(mol, atoms_data)

    atom_type_map = {atom['id']: atom['type_name'] for atom in atoms_data}
    bonds_list = []
    for bond in mol.GetBonds():
        a1_idx, a2_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        a1_id = next(item['id'] for item in atoms_data if item['idx'] == a1_idx)
        a2_id = next(item['id'] for item in atoms_data if item['idx'] == a2_idx)
        type1, type2 = atom_type_map[a1_id], atom_type_map[a2_id]
        bond_key = tuple(sorted((type1, type2)))
        order_str = BOND_ORDER_MAP.get(bond_key, "Single")
        bonds_list.append({"atoms": [a1_id, a2_id], "order": order_str})

    for atom_spec in atoms_data: del atom_spec['idx']
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
