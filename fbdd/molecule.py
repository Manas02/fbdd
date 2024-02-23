"""Molecule class."""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import BRICS

from .error import InvalidSMILESError
from .fragment import Fragment


@dataclass
class Molecule:
    smiles: str

    def __post_init__(self):
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            raise InvalidSMILESError(self.smiles)

    @property
    def mol(self) -> Chem.Mol:
        return Chem.MolFromSmiles(self.smiles)

    def fragment(self) -> list[Fragment]:
        brics_bonds: list[tuple[tuple[int, int], tuple[str, str]]] = list(BRICS.FindBRICSBonds(self.mol))
        bonds: list[tuple[tuple[int, int], tuple[str, str]]] = [(i[0], (str(n + 1), str(n + 1))) for n, i in
                                                                enumerate(brics_bonds)]
        mol_frags: Chem.Mol = BRICS.BreakBRICSBonds(self.mol, bonds=bonds)
        mol_frags: list[Chem.Mol] = Chem.GetMolFrags(mol_frags, asMols=True)
        frags: list[Fragment] = [Fragment(mol=mol) for mol in mol_frags]
        return frags

    def __hash__(self):
        return hash(self.smiles)

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
