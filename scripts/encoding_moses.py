import re
from dataclasses import dataclass

from fbdd import Molecule, Fragment


# TODO: Define behavior of this class
# TODO: Write Encoder, Decoder classes
# TODO: Fragment Augmentation
# TODO: Run over all MOSES on bcp computer
# TODO: GET STATS OVER TOKENS, make plot for frequency distribution and Data Funnel from MOSES to vocab
# TODO: Run word tokenizer over the result and train BERT over it
# TODO: Try BART, AlBERT and compare loss and performance
# TODO: Make a benchmark test for these MLMs.


class Token:
    def __init__(self, token: str):
        self.token = token
        self.valid_pattern = re.compile(r"\[\d+\*]|(\[\*])")

    def check_validity(self):
        if bool(re.search(self.valid_pattern, self.token)):
            return self.token
        return None


@dataclass
class Tokenizer:
    """
    Tokenizer takes in a Molecule object and returns a list of Tokens
    """
    mol: Molecule

    def tokenize(self) -> list[Token]:
        return NotImplemented(self.mol)


class FragmentAugmentation:
    randomize: int = 0

    def augment(self):
        if not self.randomize:
            return
        return NotImplemented()


class Encoder:
    """
    Encoder converts tokens into string with Fragment and Connectivity information separated by ` `.
    """
    tokenizer: Tokenizer
    augment: FragmentAugmentation

    def encode(self, mol: Molecule) -> str:
        return NotImplemented(self)


class AttachFragmentsAtIdentifier:
    fragment1: Fragment
    fragment2: Fragment

    def attach(self):
        return NotImplemented(self.fragment1, self.fragment2)


class Decoder:
    """
    Decoder converts string with Fragment and Connectivity information into valid SMILES.
    """
    rxn: AttachFragmentsAtIdentifier

    def decode(self, encoded_string: str) -> Molecule:
        return NotImplemented(self, encoded_string)
