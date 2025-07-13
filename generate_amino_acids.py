# generate_amino_acids.py
import json
import os
from auto_typer import build_molecule_from_smiles

# (Full Name, 3-Letter Code, 1-Letter Code, SMILES String)
AMINO_ACIDS = [
    ("Alanine", "Ala", "A", "[NH3+]C(C)C(=O)[O-]"),
    ("Glycine", "Gly", "G", "[NH3+]CC(=O)[O-]"),
    ("Isoleucine", "Ile", "I", "[NH3+]C(C(C)CC)C(=O)[O-]"),
    ("Leucine", "Leu", "L", "[NH3+]C(CC(C)C)C(=O)[O-]"),
    ("Methionine", "Met", "M", "[NH3+]C(CCSC)C(=O)[O-]"),
    ("Proline", "Pro", "P", "C1C([NH2+]C1)C(=O)[O-]"),
    ("Valine", "Val", "V", "[NH3+]C(C(C)C)C(=O)[O-]"),
    ("Asparagine", "Asn", "N", "[NH3+]C(CC(=O)N)C(=O)[O-]"),
    ("Cysteine", "Cys", "C", "[NH3+]C(CS)C(=O)[O-]"),
    ("Glutamine", "Gln", "Q", "[NH3+]C(CCC(=O)N)C(=O)[O-]"),
    ("Serine", "Ser", "S", "[NH3+]C(CO)C(=O)[O-]"),
    ("Threonine", "Thr", "T", "[NH3+]C(C(C)O)C(=O)[O-]"),
    ("Phenylalanine", "Phe", "F", "[NH3+]C(Cc1ccccc1)C(=O)[O-]"),
    ("Tryptophan", "Trp", "W", "[NH3+]C(Cc1c[nH]c2ccccc12)C(=O)[O-]"),
    ("Tyrosine", "Tyr", "Y", "[NH3+]C(Cc1ccc(O)cc1)C(=O)[O-]"),
    ("Arginine", "Arg", "R", "[NH3+]C(CCCNC(=[NH2+])N)C(=O)[O-]"),
    # explicit, non-aromatic SMILES for Histidine to prevent kekulization failure.
    ("Histidine", "His", "H", "[NH3+][C@@H](CC1=CNC=N1)C(=O)[O-]"),
    ("Lysine", "Lys", "K", "[NH3+]C(CCCC[NH3+])C(=O)[O-]"),
    ("Aspartic Acid", "Asp", "D", "[NH3+]C(CC(=O)[O-])C(=O)[O-]"),
    ("Glutamic Acid", "Glu", "E", "[NH3+]C(CCC(=O)[O-])C(=O)[O-]"),
]

def generate_and_save_molecules():
    output_dir = os.path.join("assets", "molecules")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Generating amino acid files in '{output_dir}'...")

    for name, three_letter, one_letter, smiles in AMINO_ACIDS:
        print(f"  - Processing {name} ({three_letter}/{one_letter})...")

        molecule_json = build_molecule_from_smiles(smiles, f"{name} (Zwitterion)")

        # THE FIX: Robust error handling for None return.
        if molecule_json is None:
            print(f"    [ERROR] Could not generate {name}: RDKit failed to process SMILES.")
            continue
        if "error" in molecule_json:
            print(f"    [ERROR] Could not generate {name}: {molecule_json['error']}")
            continue

        filename = f"{name.lower().replace(' ', '_')}.json"
        filepath = os.path.join(output_dir, filename)

        with open(filepath, 'w') as f:
            json.dump(molecule_json, f, indent=2)

        print(f"    [SUCCESS] Saved to {filepath}")

if __name__ == "__main__":
    generate_and_save_molecules()
    print("\nGeneration complete.")
