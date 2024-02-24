from itertools import permutations

import re
from math import factorial
import random


from rdkit import Chem
from rdkit.Chem import BRICS, AllChem


def random_enumerate(iterable):
    indices = list(range(len(iterable)))
    random.shuffle(indices)
    for i, idx in enumerate(indices):
        yield i, iterable[idx]

def to_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol)

def random_n_unique_shuffles(lst, n):
    return random.choices(list(set(permutations(lst))), k=n)

def fragment_augmentation(smi: str, augment: int = 5):
    mol: Chem.Mol = Chem.MolFromSmiles(smi)
    brics_bonds: list[tuple[tuple[int, int], tuple[str, str]]] = list(BRICS.FindBRICSBonds(mol))

    if len(brics_bonds) < 3 and augment > 2:
        augment = factorial(len(brics_bonds))
    
    if augment < 2:
        augment = 1
        
    augmented_brics_bonds = [list(i) for i in random_n_unique_shuffles([i[0] for i in brics_bonds], augment)]

    for _brics_bonds in augmented_brics_bonds:
        bonds: list[tuple[tuple[int, int], tuple[str, str]]] = [(i, (str(n + 1), str(n + 1))) 
                                                                for n, i in enumerate(_brics_bonds)]
        mol_frags: Chem.Mol = BRICS.BreakBRICSBonds(mol, bonds=bonds)
        yield to_smiles(mol_frags).split(".")


def encoder(frags: list[str]) -> str:
    if len(frags) == 1:
        raise ValueError('BRICS Failed') 
    pattern = re.compile(r'\[\d+\*\]')
    tokenized_fragment = ""
    for frag in frags:
        tokenized_fragment += re.sub(pattern, '[*]', frag) + " " + " ".join(pattern.findall(frag)) + " "
    return tokenized_fragment[:-1]

def extract_attachment_token(fragment):
    return re.findall(r'\d+\*', fragment)

def extract_fragments(encoded_fragment) -> str:
    return ' '.join(item for item in encoded_fragment.split() if not re.compile(r'\[\d+\*\]').match(item))


def combine_fragments(smi1:str, smi2:str, attachment_token:str) -> str|ValueError:
    # Fails to combine fragments with double bond
    m1 = Chem.MolFromSmiles(smi1) 
    m2 = Chem.MolFromSmiles(smi2) 
    rxn = AllChem.ReactionFromSmarts(f"[{attachment_token}][*:1].[{attachment_token}][*:2]>>[*:1][*:2]")
    results = rxn.RunReactants([m1, m2])
    if results:
        return Chem.MolToSmiles(results[0][0])
    raise ValueError(smi1, smi2)
    

def _decode(encoded_fragment:str) -> str | ValueError:
        s: str = extract_fragments(encoded_fragment)
        replacements: list = extract_attachment_token(encoded_fragment)
        if not replacements: 
            raise ValueError('Molecule was not fragmented')
        
        # Calculate the length of the replacement string
        replacement_length = len(replacements[0])
        index = 0
    
        # Iterate over the replacements and replace each occurrence sequentially
        while '[*]' in s:
            # Find the index of the next occurrence of '[*]'
            start_index = s.find('[*]', index)
            if start_index == -1:
                break  # Exit the loop if no more occurrences are found
    
            # Determine the replacement value
            replacement = replacements.pop(0)
    
            # Replace '[*]' with the replacement value
            s = s[:start_index+1] + replacement + s[start_index + replacement_length:]
    
            # Update the index to continue searching for the next occurrence
            index = start_index + len(replacement)
    
        return s   


def decoder(encoded_fragment):
        frags = _decode(encoded_fragment).split()
        base_fragment = frags[0]
        frags = frags[1:]
        
        condition = True
        while condition:
            for attachment_token in extract_attachment_token(base_fragment):
                for n, frag in enumerate(frags):
                    if attachment_token in frag:
                        base_fragment = combine_fragments(base_fragment, frag, attachment_token)
                        frags.remove(frag)
                        break  # Move to the next attachment token once one is found and processed
                        
            if not frags:
                break
            
            if not extract_attachment_token(base_fragment):
                break
            
        return base_fragment