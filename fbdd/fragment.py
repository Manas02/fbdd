"""Fragment class."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from .error import InvalidSMARTSError


@dataclass
class Fragment:
    mol: Chem.Mol

    def to_smarts(self) -> str | InvalidSMARTSError:
        smarts = Chem.MolToSmarts(self.mol)
        if smarts:
            return smarts
        raise InvalidSMARTSError()
