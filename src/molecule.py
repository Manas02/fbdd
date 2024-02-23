"""Molecule class definition."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import BRICS

from error import InvalidSMILESError
from fragment import Fragment, Fragments


@dataclass
class Molecule:
    mol: Chem.Mol

    @staticmethod
    def to_smiles(mol) -> str | InvalidSMILESError:
        smiles = Chem.MolToSmiles(mol)
        if smiles:
            return smiles
        raise InvalidSMILESError(mol)

    def smiles(self) -> str:
        return self.to_smiles(self.mol)

    def fragment(self) -> Fragments:
        brics_bonds: list[tuple[tuple[int, int], tuple[str, str]]] = list(BRICS.FindBRICSBonds(self.mol))
        bonds: list[tuple[tuple[int, int], tuple[str, str]]] = [(i[0], (str(n + 1), str(n + 1))) for n, i in
                                                                enumerate(brics_bonds)]
        mol_frags: Chem.Mol = BRICS.BreakBRICSBonds(self.mol, bonds=bonds)
        mol_frags: list[Chem.Mol] = Chem.GetMolFrags(mol_frags, asMols=True)
        frags: Fragments = [Fragment(mol=mol) for mol in mol_frags]
        return frags

    def __hash__(self):
        return hash(self.to_smiles(self.mol))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()


Molecules = list[Molecule]
