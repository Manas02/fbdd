"""Fragment class."""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from .error import InvalidSMARTSError


@dataclass
class Fragment:
    mol: Chem.Mol

    def to_smarts(self) -> str | InvalidSMARTSError:
        try:
            smarts = Chem.MolToSmarts(self.mol)
            return smarts
        except Exception as e:
            raise InvalidSMARTSError(e)
