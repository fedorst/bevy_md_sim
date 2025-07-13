# test_autotyper.py
import unittest
import json
from collections import Counter
import os

# Import the function we want to test from your script
from auto_typer import build_molecule_from_smiles

# --- NEW: Load the force field once to check against it ---
# This assumes the script is run from the project's root directory.
FORCE_FIELD_PATH = os.path.join('assets', 'force_field.json')
with open(FORCE_FIELD_PATH) as f:
    FORCE_FIELD = json.load(f)

# Create a set of valid bonds for fast lookups.
# We store them as a canonical tuple: (sorted_type_1, sorted_type_2, order)
VALID_BONDS = set()
for bond_data in FORCE_FIELD['bonds']:
    types = tuple(sorted(bond_data['types']))
    # Default to "Single" if 'order' key is missing, for backwards compatibility
    order = bond_data.get('order', 'Single')
    VALID_BONDS.add((types[0], types[1], order))

class TestAtomTyper(unittest.TestCase):
    """
    A suite of unit tests to verify the correctness of the atom typing logic.
    Each test provides a SMILES string for a known molecule and asserts that
    the generated atom types match our expectations.
    """

    # --- NEW: A more comprehensive validation helper ---
    def _validate_molecule(self, smiles, name):
        """
        Runs the molecule builder and performs critical validation checks:
        1. Fails if any atom is not assigned a specific type (e.g., remains 'C' or 'H').
        2. Fails if any generated bond is not defined in the force field.
        """
        mol_json = build_molecule_from_smiles(smiles, name)
        self.assertIsNotNone(mol_json, f"Molecule generation failed for {name} ({smiles})")
        self.assertNotIn("error", mol_json, f"SMILES parsing returned an error for {name}: {mol_json.get('error')}")

        atoms = mol_json['atoms']
        bonds = mol_json['bonds']

        # 1. Check for any untyped atoms
        untyped_atoms = [atom for atom in atoms if atom['type_name'] in ['C', 'H', 'O', 'N']]
        self.assertEqual(len(untyped_atoms), 0, f"Found untyped atoms for {name}: {untyped_atoms}")

        # 2. Check for any undefined bonds
        atom_type_map = {atom['id']: atom['type_name'] for atom in atoms}
        undefined_bonds = []
        for bond in bonds:
            atom1_id, atom2_id = bond['atoms']
            order = bond.get('order', 'Single')

            type1 = atom_type_map[atom1_id]
            type2 = atom_type_map[atom2_id]

            bond_key = tuple(sorted((type1, type2)))
            full_bond_key = (bond_key[0], bond_key[1], order)

            if full_bond_key not in VALID_BONDS:
                undefined_bonds.append(full_bond_key)

        self.assertEqual(len(undefined_bonds), 0, f"Found undefined bonds for {name}: {undefined_bonds}")

        return mol_json

    def _run_test_and_get_type_counts(self, smiles, name):
        """Helper function to run the molecule builder and return type counts."""
        result_json = build_molecule_from_smiles(smiles, name)
        self.assertIsNotNone(result_json, f"Failed to generate molecule for SMILES: {smiles}")
        type_names = [atom['type_name'] for atom in result_json['atoms']]
        actual_types = Counter(type_names)

        self.addTypeEqualityFunc(dict, lambda a, b, msg=None: self.assertDictEqual(a, b, msg=f"Expected != Actual.\nActual was: {b}"))
        return actual_types

    # --- NEW: Test cases for the amino acids we'll use in the peptide builder ---

    def test_alanine_zwitterion_validation(self):
        """Validates that all atoms and bonds in Alanine are correctly defined."""
        smiles = "[NH3+]C(C)C(=O)[O-]"
        self._validate_molecule(smiles, "Alanine (Zwitterion)")

    def test_serine_zwitterion_validation(self):
        """Validates that all atoms and bonds in Serine are correctly defined."""
        smiles = "[NH3+]C(CO)C(=O)[O-]"
        self._validate_molecule(smiles, "Serine (Zwitterion)")

    def test_valine_zwitterion_validation(self):
        """Validates that all atoms and bonds in Valine are correctly defined."""
        smiles = "[NH3+]C(C(C)C)C(=O)[O-]"
        self._validate_molecule(smiles, "Valine (Zwitterion)")

    def test_aspartic_acid_zwitterion_validation(self):
        """Validates that Aspartic Acid (with its own carboxyl group) is correct."""
        smiles = "[NH3+]C(CC(=O)O)C(=O)[O-]"
        self._validate_molecule(smiles, "Aspartic Acid (Zwitterion)")

    # --- Existing tests below (can be kept for regression) ---

    def test_water(self):
        """Tests a simple water molecule."""
        self._validate_molecule("O", "Water")

    def test_ethanol(self):
        """Tests ethanol, which has two different carbons and a hydroxyl group."""
        self._validate_molecule("CCO", "Ethanol")

    def test_acetaldehyde(self):
        """Tests acetaldehyde, featuring a carbonyl group."""
        self._validate_molecule("CC=O", "Acetaldehyde")

    def test_methylamine(self):
        self._validate_molecule("CN", "Methylamine")

    def test_glycine_zwitterion(self):
        # The validation test for Glycine is implicitly handled by the old test.
        # We can use the more powerful helper here too.
        smiles = "C(C(=O)[O-])[NH3+]"
        self._validate_molecule(smiles, "Glycine (Zwitterion)")
        # You can optionally keep the count check if you want to be extra sure.
        expected_types = { "CT": 1, "C_CO": 1, "O_COO": 2, "N_PLUS": 1, "H_C": 2, "H_N": 3 }
        mol_json = build_molecule_from_smiles(smiles, "Glycine (Zwitterion)")
        type_names = [atom['type_name'] for atom in mol_json['atoms']]
        self.assertDictEqual(expected_types, Counter(type_names))

    def test_phenylalanine_zwitterion(self):
        smiles = "c1ccc(C[C@@H]([NH3+])C(=O)[O-])cc1"
        self._validate_molecule(smiles, "Phenylalanine (Zwitterion)")

    def test_arginine_zwitterion(self):
        smiles = "[NH3+][C@@H](CCCNC(=N)N)C(=O)[O-]"
        self._validate_molecule(smiles, "Arginine (Zwitterion)")

    def test_nitromethane(self):
        smiles = "C[N+](=O)[O-]"
        self._validate_molecule(smiles, "Nitromethane")


if __name__ == '__main__':
    unittest.main()
