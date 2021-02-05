# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

from typing import Sequence

from rdkit.Chem.rdchem import Mol


def clear_atom_mapping(mols: Sequence[Mol]) -> None:
    """
    Remove the mapping information in RDKit Mol objects.
    """
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                atom.ClearProp('molAtomMapNumber')
