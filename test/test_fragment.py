from fbdd import Fragment, InvalidSMARTSError
from rdkit import Chem
import pytest


@pytest.fixture
def test_fragment():
    mol = Chem.MolFromSmiles('C')
    fragment = Fragment(mol)
    return fragment


def test_to_smarts_valid(test_fragment):
    smarts = test_fragment.to_smarts()
    assert smarts == '[#6]'


def test_to_smarts_invalid():
    mol = Chem.MolFromSmiles('')
    fragment = Fragment(mol)
    with pytest.raises(InvalidSMARTSError):
        fragment.to_smarts()
