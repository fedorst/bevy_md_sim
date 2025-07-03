# auto_typer.py (Version 2)
import json
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem


# Keep the same rules as before
TYPE_RULES = [
    ("N_PLUS", "[N+1;H3;D3;X3]"), ("O_COO", "[O-1;D1;X1]"), ("C_GUA", "[C;X3](=[N;X2])([N;X3])[N;X3]"),
    ("N_GUA", "[N;X3]([C;X3](=[N;X2]))"), ("CA", "c"), ("HA", "[H]c"), ("C_CO", "[C;X3]=[O;X1]"),
    ("O_CO", "[O;X1]=[C;X3]"), ("N_AMINE", "[N;H2;D3;X3]"), ("O_H", "[O;H1;D2;X2]"),
    ("H_O", "[H][O;H1;D2]"), ("CT", "[C;X4]"), ("H_C", "[H][C;X4]"),
]

def build_molecule_from_smiles(smiles_string, name):
    """Generates a full molecule JSON from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol:
        return None

    # Add hydrogens, generate 3D coords, and do a quick optimization
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    # --- Start building the JSON structure ---
    atoms_data = []
    # Use the conformer to get 3D positions
    conformer = mol.GetConformer()
    for atom in mol.GetAtoms():
        pos = conformer.GetAtomPosition(atom.GetIdx())
        # Generate a unique ID like C1, C2, H3, etc.
        atom_id = f"{atom.GetSymbol()}{atom.GetIdx() + 1}"
        atoms_data.append({
            "id": atom_id,
            "type_name": "", # This will be filled in next
            "element": atom.GetSymbol(),
            "pos": [pos.x / 10.0, pos.y / 10.0, pos.z / 10.0] # Convert Angstrom to nm
        })

    # --- Assign types using SMARTS rules ---
    assigned_types = {}


    # Nitro group special case
    nitro_pat = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pat)
    for match in nitro_matches:
        n_idx, o1_idx, o2_idx = match
        assigned_types[n_idx] = "N_NO2"
        assigned_types[o1_idx] = "O_NO2"
        assigned_types[o2_idx] = "O_NO2"

    for type_name, smarts in TYPE_RULES:
        pattern = Chem.MolFromSmarts(smarts)
        if not pattern: continue
        matches = mol.GetSubstructMatches(pattern)
        for match_indices in matches:
            for idx in match_indices:
                if idx not in assigned_types:
                    assigned_types[idx] = type_name

    for i, atom_spec in enumerate(atoms_data):
        if i in assigned_types:
            atom_spec['type_name'] = assigned_types[i]
        else:
            atom_spec['type_name'] = atom_spec['element']


    # --- Extract bonds ---
    bonds_data = []
    double_bonds_data = []
    for bond in mol.GetBonds():
        a1_id = f"{mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()}{bond.GetBeginAtomIdx() + 1}"
        a2_id = f"{mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()}{bond.GetEndAtomIdx() + 1}"

        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonds_data.append([a1_id, a2_id])
        else: # Treat single, triple, aromatic as single for this sim
            bonds_data.append([a1_id, a2_id])

    return {
        "name": name,
        "atoms": atoms_data,
        "bonds": bonds_data,
        "double_bonds": double_bonds_data
    }


def main():
    parser = argparse.ArgumentParser(description="Generate a typed molecule JSON from a SMILES string.")
    parser.add_argument("--smiles", required=True, help="SMILES string of the molecule to generate.")
    parser.add_argument("--name", default="Generated Molecule", help="Name for the molecule.")
    args = parser.parse_args()

    molecule_json = build_molecule_from_smiles(args.smiles, args.name)

    if molecule_json:
        # Print the final JSON to standard output
        print(json.dumps(molecule_json, indent=2))
    else:
        # Print an error to stderr if something goes wrong
        import sys
        print(json.dumps({"error": f"Could not parse SMILES: {args.smiles}"}), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
