import sys
sys.path.append("../src")

import unittest
from unittest.mock import MagicMock
from rdkit import Chem

from molecule import Molecule
from error import InvalidSMILESError
from fragment import Fragment, Fragments


class TestMolecule(unittest.TestCase):

    def setUp(self):
        # Create a test molecule
        self.test_smiles = 'CCO'  # Ethanol
        self.test_mol = Chem.MolFromSmiles(self.test_smiles)
        self.molecule = Molecule(self.test_mol)

    def test_to_smiles_valid(self):
        # Test if the to_smiles method returns valid SMILES
        self.assertEqual(self.molecule.to_smiles(self.test_mol), self.test_smiles)

    def test_to_smiles_invalid(self): # FIXME: here
        # Test if the to_smiles method raises InvalidSMILESError for invalid molecule
        invalid_mol = MagicMock()
        invalid_mol.return_value = None  # Mocking a molecule without valid SMILES
        with self.assertRaises(InvalidSMILESError):
            self.molecule.to_smiles(invalid_mol)

    def test_fragment(self):
        # Test the fragment method
        fragments = self.molecule.fragment()
        self.assertIsInstance(fragments, Fragments) # FIXME: due to UnFragmentableMolecule
        self.assertTrue(all(isinstance(fragment, Fragment) for fragment in fragments))

    def test_hash_equality(self):
        # Test hash and equality methods
        mol_1 = Molecule(Chem.MolFromSmiles('CCO'))
        mol_2 = Molecule(Chem.MolFromSmiles('CCO'))
        mol_3 = Molecule(Chem.MolFromSmiles('CCC'))

        self.assertEqual(hash(mol_1), hash(mol_2))
        self.assertNotEqual(hash(mol_1), hash(mol_3))
        self.assertEqual(mol_1, mol_2)
        self.assertNotEqual(mol_1, mol_3)


if __name__ == '__main__':
    unittest.main()
