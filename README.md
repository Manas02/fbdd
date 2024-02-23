# FragmentBERT

## Molecule and Fragment Classes

This repository contains Python code for the `Molecule` and `Fragment` classes, which are designed to handle molecules and their fragments in the context of computational chemistry.

### Molecule Class

The `Molecule` class represents a chemical molecule and provides functionality for handling molecule-related operations such as converting molecules to SMILES (Simplified Molecular Input Line Entry System) notation, fragmenting molecules, and hashing/equality comparison.

### Features

- Convert a molecule to its SMILES representation.
- Fragment a molecule into smaller molecular fragments.
- Support for hashing and equality comparison based on molecular structure.

### Usage

```python
from fbdd import Molecule

# Create a molecule from a SMILES string
smi = 'CCO'
molecule = Molecule(smi)

# Get the SMILES representation of the molecule
mol = molecule.mol
print("Molecule:", mol)

# Fragment the molecule into smaller fragments
fragments = molecule.fragment()
for fragment in fragments:
    print("Fragment SMARTS:", fragment.to_smarts())
```

### Fragment Class

The `Fragment` class represents a fragment of a chemical molecule and provides functionality for handling fragment-related operations such as converting fragments to SMARTS notation.

### Features

- Convert a fragment to its SMARTS representation.

### Usage

```python
from rdkit import Chem
from fbdd import Fragment

# Create a fragment from a SMILES string
mol = Chem.MolFromSmiles('C')
fragment = Fragment(mol)

# Get the SMARTS representation of the fragment
smarts = fragment.to_smarts()
print("SMARTS:", smarts)
```
