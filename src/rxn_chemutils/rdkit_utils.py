# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

from typing import Iterable

from rdkit.Chem import CombineMols
from rdkit.Chem.rdchem import Mol


def clear_atom_mapping(mols: Iterable[Mol]) -> None:
    """
    Remove the mapping information in RDKit Mol objects.
    """
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                atom.ClearProp("molAtomMapNumber")


def combine_mols(mols: Iterable[Mol]) -> Mol:
    """Combine multiple RDKit Mols into one.

    RDKit has a function for that but it only supports two mols.

    Args:
        mols: RDKit Mols to combine.
    """

    resulting_mol = Mol()
    for mol in mols:
        resulting_mol = CombineMols(resulting_mol, mol)
    return resulting_mol
