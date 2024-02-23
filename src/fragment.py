"""Fragment class definition."""
from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from error import InvalidSMARTSError


@dataclass
class Fragment:
    mol: Chem.Mol

    def to_smarts(self) -> str|InvalidSMARTSError:
        return Chem.MolToSmarts(self.mol)


Fragments = list[Fragment]