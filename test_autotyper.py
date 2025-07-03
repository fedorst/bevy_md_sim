# test_autotyper.py
import unittest
import json
from collections import Counter

# Import the function we want to test from your script
from auto_typer import build_molecule_from_smiles

class TestAtomTyper(unittest.TestCase):
    """
    A suite of unit tests to verify the correctness of the atom typing logic.
    Each test provides a SMILES string for a known molecule and asserts that
    the generated atom types match our expectations.
    """


    def _run_test_and_get_type_counts(self, smiles, name):
        """Helper function to run the molecule builder and return type counts."""
        result_json = build_molecule_from_smiles(smiles, name)
        self.assertIsNotNone(result_json, f"Failed to generate molecule for SMILES: {smiles}")
        type_names = [atom['type_name'] for atom in result_json['atoms']]
        actual_types = Counter(type_names)

        # --- LOGGING on Failure ---
        # Add a custom failure message that includes the actual result.
        # This makes debugging much faster.
        self.addTypeEqualityFunc(dict, lambda a, b, msg=None: self.assertDictEqual(a, b, msg=f"Expected != Actual.\nActual was: {b}"))

        return actual_types

    def test_water(self):
        """Tests a simple water molecule."""
        smiles = "O"
        expected_types = {
            "O_H": 1,  # The oxygen in water
            "H_O": 2   # The two hydrogens on water
        }
        actual_types = self._run_test_and_get_type_counts(smiles, "Water")
        self.assertDictEqual(expected_types, actual_types)

    def test_ethanol(self):
        """Tests ethanol, which has two different carbons and a hydroxyl group."""
        smiles = "CCO"
        expected_types = {
            "CT": 2,   # Two sp3 carbons
            "H_C": 5,  # Five hydrogens on carbons
            "O_H": 1,  # One hydroxyl oxygen
            "H_O": 1   # One hydrogen on the hydroxyl oxygen
        }
        actual_types = self._run_test_and_get_type_counts(smiles, "Ethanol")
        self.assertDictEqual(expected_types, actual_types)

    def test_acetaldehyde(self):
        """Tests acetaldehyde, featuring a carbonyl group."""
        smiles = "CC=O"
        expected_types = {
            "CT": 1,   # The methyl carbon
            "C_CO": 1, # The carbonyl carbon
            "O_CO": 1, # The carbonyl oxygen
            "H_C": 4   # 3 on the methyl group, 1 on the aldehyde group
        }
        actual_types = self._run_test_and_get_type_counts(smiles, "Acetaldehyde")
        self.assertDictEqual(expected_types, actual_types)

    def test_methylamine(self):
        smiles = "CN"
        expected_types = {
            "CT": 1,
            "N_AMINE": 1,
            "H_C": 3,
            "H_N": 2,  # <-- CORRECTED: Was H_O
        }
        actual_types = self._run_test_and_get_type_counts(smiles, "Methylamine")
        self.assertDictEqual(expected_types, actual_types)

    def test_glycine_zwitterion(self):
        smiles = "C(C(=O)[O-])[NH3+]"
        expected_types = {
            "CT": 1,
            "C_CO": 1,
            "O_COO": 2,
            "N_PLUS": 1,
            "H_C": 2,
            "H_N": 3,  # <-- CORRECTED: Was H_O
        }
        actual_types = self._run_test_and_get_type_counts(smiles, "Glycine (Zwitterion)")
        self.assertDictEqual(expected_types, actual_types)

    def test_phenylalanine_zwitterion(self):
        smiles = "c1ccc(C[C@@H]([NH3+])C(=O)[O-])cc1"
        expected_types = {
            "CA": 6,
            "HA": 5,
            "CT": 2,
            "C_CO": 1,
            "O_COO": 2,
            "N_PLUS": 1,
            "H_C": 3,
            "H_N": 3,  # <-- CORRECTED: Was H_O
        }
        actual_types = self._run_test_and_get_type_counts(smiles, "Phenylalanine (Zwitterion)")
        self.assertDictEqual(expected_types, actual_types)

    def test_arginine_zwitterion(self):
        """
        Tests arginine, which has a complex guanidinium group.
        We test the zwitterionic form.
        """
        smiles = "[NH3+][C@@H](CCCNC(=N)N)C(=O)[O-]"
        expected_types = {
            # Backbone
            "N_PLUS": 1,
            "CT": 1,      # Alpha-carbon
            "C_CO": 1,
            "O_COO": 2,
            "H_N": 3,     # On the N+
            "H_C": 1,     # On the alpha-carbon
            # Side Chain
            "CT": 3,      # Beta, Gamma, Delta carbons
            "H_C": 6,     # 2 on each of the 3 CTs
            "N_GUA": 3,   # 3 nitrogens in the guanidinium group
            "C_GUA": 1,   # 1 carbon in the guanidinium group
            "H_N": 4      # 4 hydrogens on the guanidinium nitrogens
        }

        # We need to combine the counts
        final_expected = {
            "N_PLUS": 1, "CT": 4, "C_CO": 1, "O_COO": 2,
            "C_GUA": 1, "N_GUA": 3,
            "H_N": 7, "H_C": 7,
        }

        actual_types = self._run_test_and_get_type_counts(smiles, "Arginine (Zwitterion)")
        self.assertDictEqual(final_expected, actual_types)

    def test_nitromethane(self):
        """Tests a nitro group, which is a special-cased functional group."""
        smiles = "C[N+](=O)[O-]"
        # You would need to add N_NO2 and O_NO2 to your force field for this to be useful
        expected_types = {
            "CT": 1,
            "H_C": 3,
            "N_NO2": 1,
            "O_NO2": 2
        }
        actual_types = self._run_test_and_get_type_counts(smiles, "Nitromethane")
        self.assertDictEqual(expected_types, actual_types)

if __name__ == '__main__':
    unittest.main()
