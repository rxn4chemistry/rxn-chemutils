import typing
from collections import Counter
from typing import List

from rdkit.Chem import MolFromSmiles, Mol, Atom, AddHs

from rxn_chemutils.conversion import InvalidSmiles, canonicalize_smiles


def is_valid_smiles(s: str) -> bool:
    """
    Whether a given string corresponds to a valid SMILES string.

    Args:
        s: string to check.

    Returns:
        True if s is a valid SMILES string, else False.
    """
    try:
        canonicalize_smiles(s)
        return True
    except InvalidSmiles:
        return False


def equivalent_smiles(*smiles: str) -> bool:
    """
    Returns true if all the given SMILES strings are equivalent.
    Will catch the exceptions for invalid SMILES and return false in that case.
    """
    try:
        canonical_smiles = [canonicalize_smiles(s) for s in smiles]
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
