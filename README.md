# RXN chemistry utilities package

[![Actions tests](https://github.com/rxn4chemistry/rxn-chemutils/actions/workflows/tests.yaml/badge.svg)](https://github.com/rxn4chemistry/rxn-chemutils/actions)

This repository contains various chemistry-related Python utilities used in the RXN universe.
For general utilities not related to chemistry, see our other repository [`rxn-utilities`](https://github.com/rxn4chemistry/rxn-utilities).

Links:
* [GitHub repository](https://github.com/rxn4chemistry/rxn-chemutils)
* [Documentation](https://rxn4chemistry.github.io/rxn-chemutils/)
* [PyPI package](https://pypi.org/project/rxn-chem-utils/)

## System Requirements

This package is supported on all operating systems. 
It has been tested on the following systems:
+ macOS: Big Sur (11.1)
+ Linux: Ubuntu 18.04.4

A Python version of 3.7 or greater is recommended.

## Installation guide

The package can be installed from Pypi:
```bash
pip install rxn-chem-utils
```

For local development, the package can be installed with:
```bash
pip install -e .[dev]
```

The `RDKit` dependency is not installed automatically and can be installed via Conda or Pypi:
```bash
# Install RDKit from Conda
conda install -c conda-forge rdkit

# Install RDKit from Pypi
pip install rdkit
```

## Package highlights

### Convert between compound representations

There are functions to convert between SMILES, `RDKit.Mol`, MDL, InChI, etc. 
All of them work in a similar way:
```pycon
>>> from rxn.chemutils.conversion import smiles_to_mol, mol_to_smiles
>>> mol = smiles_to_mol("CO(C)")
>>> mol_to_smiles(mol)
'COC'
```

The functions raise exceptions when failing, and allow to be used without sanitization.
```pycon
>>> mol = smiles_to_mol("CFC")
Traceback (most recent call last):
[...]
rxn.chemutils.exceptions.InvalidSmiles: "CFC" is not a valid SMILES string
```
```pycon
>>> mol = smiles_to_mol("CFC", sanitize=False)
>>> mol_to_smiles(mol)
'CFC'
```

### Reaction SMILES

The package supports different kinds of reaction SMILES, which, internally, are stored as [`ReactionEquation`s](./src/rxn/chemutils/reaction_equation.py).

To convert to and from `ReactionEquation`, a few functions are provided:
* [`parse_reaction_smiles`]() and [`to_reaction_smiles`](), if you know already the format.
* [`parse_any_reaction_smiles`](), if you don't know the format or want to be flexible.

Examples:
```pycon
>>> from rxn.chemutils.reaction_smiles import ReactionFormat, determine_format, parse_reaction_smiles, to_reaction_smiles, parse_any_reaction_smiles
>>> rxn_smiles = "CC.O.[Na+]~[Cl-]>>CCO"
>>> determine_format(rxn_smiles)
<ReactionFormat.STANDARD_WITH_TILDE: 3>
>>> parse_reaction_smiles(rxn_smiles, ReactionFormat.STANDARD_WITH_TILDE)
ReactionEquation(reactants=['CC', 'O', '[Na+].[Cl-]'], agents=[], products=['CCO'])
>>> parse_any_reaction_smiles(rxn_smiles)
ReactionEquation(reactants=['CC', 'O', '[Na+].[Cl-]'], agents=[], products=['CCO'])
>>> to_reaction_smiles(parse_any_reaction_smiles(rxn_smiles), ReactionFormat.EXTENDED)
'CC.O.[Na+].[Cl-]>>CCO |f:2.3|'
```

### Multicomponent SMILES

Sometimes, it is necessary to represent multiple compounds as one single SMILES string.
For fragments / ions, it becomes necessary to distinguish between what parts belong together as one compound, and what are differentt compounds.
In such "multitcomponent SMILES", we typically use tildes, `~`, to indicate that different SMILES fragments belong to the same compound.

```pycon
>>> from rxn.chemutils.multicomponent_smiles import multicomponent_smiles_to_list, list_to_multicomponent_smiles
>>> list_to_multicomponent_smiles(["CC", "[Na+].[Cl-"], fragment_bond="~")
'CC.[Na+]~[Cl-'
>>> multicomponent_smiles_to_list('CC.[Na+]~[Cl-', fragment_bond="~")
['CC', '[Na+].[Cl-']
```

### Canonicalization

Canonicalization of compounds, with the possibility to remove the valence check:
```pycon
>>> from rxn.chemutils.conversion import canonicalize_smiles
>>> canonicalize_smiles("CC(O)")
'CCO'
>>> canonicalize_smiles("ABCD")  # Invalid SMILES
Traceback (most recent call last):
[...]
rxn.chemutils.exceptions.InvalidSmiles: "ABCD" is not a valid SMILES string
>>> canonicalize_smiles("CF(C)")  # Invalid valence, fails by default
Traceback (most recent call last):
[...]
rxn.chemutils.exceptions.InvalidSmiles: "CFC" is not a valid SMILES string
>>> canonicalize_smiles("CF(C)", check_valence=False)  # Invalid valence, does not fail
'CFC'
```

Canonicalization of any kind of SMILES (components, multicomponent SMILES, reaction SMILES, etc.), again with the possibility to disable the valence check.
Note that the resulting string is in the same format.
```pycon
>>> from rxn.chemutils.miscellaneous import canonicalize_any
>>> canonicalize_any()
>>> canonicalize_any("[Na+].[Cl-]")
'[Cl-].[Na+]'
>>> canonicalize_any("OC.C(O)~CF(C)", check_valence=False)
'CO.CFC~CO'
>>> canonicalize_any("CC(C)>C(O)>C(O)")
'CCC>CO>CO'
>>> canonicalize_any("CO.O.C>>C(O) |f:1.2|")
'CO.C.O>>CO |f:1.2|'
```

The executable `rxn-canonicalize` (installed with the package), which works either on files or on stdin
```shell
rxn-canonicalize --help
```

### Augmentation

See [`smiles_randomization.py`](./src/rxn/chemutils/smiles_randomization.py) and [`smiles_augmenter.py`](./src/rxn/chemutils/smiles_augmenter.py) for the augmentation of compound SMILES and reaction SMILES strings.

### Others

Without going into details, the package also does the following:
* Tokenization and detokenization of SMILES strings in [`tokenization.py`](./src/rxn/chemutils/tokenization.py), and the executables `rxn-tokenize` and `rxn-detokenize`.
* Easy combination of precursor SMILES and product SMILES into a reaction SMILES with the [`ReactionCombiner`](./src/rxn/chemutils/reaction_combiner.py), and the executable `rxn-combine-reaction`.
* Parsing of RDFs into reaction SMILES: different [modules](./src/rxn/chemutils/rdf), and the executable `rxn-rdf-to-smiles`.
* ... and many others.