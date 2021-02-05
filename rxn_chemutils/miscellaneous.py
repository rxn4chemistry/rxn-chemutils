# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

import re
import typing
from collections import Counter
from typing import List

from rdkit.Chem import MolFromSmiles, Mol, Atom, AddHs

from rxn_chemutils.conversion import canonicalize_smiles
from rxn_chemutils.exceptions import InvalidSmiles


def is_valid_smiles(smiles: str, check_valence: bool = True) -> bool:
    """
    Whether a given string corresponds to a valid SMILES string.

    Args:
        smiles: string to check.
        check_valence: whether to check the valence.

    Returns:
        True if the given SMILES is valid, else False.
    """
    try:
        canonicalize_smiles(smiles, check_valence=check_valence)
        return True
    except InvalidSmiles:
        return False


def equivalent_smiles(*smiles: str, check_valence: bool = False) -> bool:
    """
    Returns true if all the given SMILES strings are equivalent.
    Will catch the exceptions for invalid SMILES and return false in that case.

    Args:
        smiles: multiple SMILES to check for equivalence.
        check_valence: if True, molecules with invalid valence will be invalidated.
    """
    try:
        canonical_smiles = [canonicalize_smiles(s, check_valence=check_valence) for s in smiles]
        return len(set(canonical_smiles)) == 1
    except InvalidSmiles:
        return False


def atom_type_counter(smiles: str) -> typing.Counter[str]:
    """
    Return a counter of atom types (as symbols).
    """

    mol: Mol = AddHs(MolFromSmiles(smiles))
    atoms: List[Atom] = mol.GetAtoms()
    return Counter(atom.GetSymbol() for atom in atoms)


def remove_atom_mapping(smiles: str) -> str:
    """
    Remove the atom mapping of a reaction SMILES.

    The resulting SMILES strings will still contain brackets and it may be
    advisable to canonicalize them as a postprocessing step.

    Args:
        smiles: SMILES string potentially containing mapping information.

    Returns:
        A SMILES string without atom mapping information.
    """

    # We look for ":" followed by digits before a "]"
    return re.sub(r':\d+]', ']', smiles)
