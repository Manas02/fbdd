import pytest
from fbdd import Molecule, Fragment, InvalidSMILESError


@pytest.fixture
def test_molecule():
    smiles = 'CCO'
    return Molecule(smiles)


def test_fragment(test_molecule):
    fragments = test_molecule.fragment()
    assert all(isinstance(fragment, Fragment) for fragment in fragments)


def test_to_smiles_valid(test_molecule):
    assert test_molecule.smiles == 'CCO'


def test_to_smiles_invalid():
    with pytest.raises(InvalidSMILESError):
        Molecule('invalid_smiles')


def test_hash_equality():
    mol_1 = Molecule('CCO')
    mol_2 = Molecule('CCO')
    mol_3 = Molecule('CCC')

    assert hash(mol_1) == hash(mol_2)
    assert hash(mol_1) != hash(mol_3)
    assert mol_1 == mol_2
    assert mol_1 != mol_3
