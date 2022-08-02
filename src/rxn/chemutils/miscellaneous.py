import re
import typing
from collections import Counter
from functools import partial
from typing import Callable, List

from rdkit.Chem import AddHs, Atom, Mol

from .conversion import canonicalize_smiles, smiles_to_mol
from .exceptions import InvalidSmiles
from .multicomponent_smiles import (
    apply_to_multicomponent_smiles,
    sort_multicomponent_smiles,
)
from .reaction_equation import apply_to_compounds, sort_compounds
from .reaction_smiles import determine_format, parse_reaction_smiles, to_reaction_smiles

CHIRAL_CENTER_PATTERN = re.compile(
    r"\[([^],@]+)@+([^]]*)]"
)  # Matches stereo centres, and groups what comes before and after "@"


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
        canonical_smiles = [
            canonicalize_smiles(s, check_valence=check_valence) for s in smiles
        ]
        return len(set(canonical_smiles)) == 1
    except InvalidSmiles:
        return False


def atom_type_counter(smiles: str) -> typing.Counter[str]:
    """
    Return a counter of atom types (as symbols).
    """

    mol: Mol = AddHs(smiles_to_mol(smiles, sanitize=False))
    atoms: List[Atom] = mol.GetAtoms()
    return Counter(atom.GetSymbol() for atom in atoms)


def remove_chiral_centers(smiles: str) -> str:
    """
    Return SMILES where all the chiral centres are removed.

    Args:
        smiles: non-atom-mapped SMILES string.

    Returns:
        SMILES with no chiral information. It is not canonical.
    """
    return re.sub(CHIRAL_CENTER_PATTERN, r"[\g<1>\g<2>]", smiles)


def remove_double_bond_stereochemistry(smiles: str) -> str:
    """
    Return SMILES where all the E/Z information on double bonds is removed.

    Args:
        smiles: SMILES string.

    Returns:
        SMILES with no sterochemical information for double bonds. The SMILES
        is not guaranteed to be canonical.
    """
    return smiles.replace("/", "").replace("\\", "")


def apply_to_any_smiles(any_smiles: str, fn: Callable[[str], str]) -> str:
    """
    Apply a given function to individual compound SMILES strings given in any kind
    of SMILES string (molecule SMILES, multicomponent SMILES, reaction SMILES).

    In the case of reaction SMILES, the format is kept.

    Args:
        any_smiles: any kind of SMILES string.
        fn: callback to apply to every compound SMILES.

    Raises:
        Exception: different kinds of exception may be raised during parsing,
            or during execution of the callback.

    Returns:
        the new (molecule, multicomponent, or reaction) SMILES string after
        application of the callback to all the component SMILES.
    """
    if ">" in any_smiles:
        # we have a reaction SMILES
        reaction_format = determine_format(any_smiles)
        reaction = parse_reaction_smiles(any_smiles, reaction_format)
        reaction = apply_to_compounds(reaction, fn)
        return to_reaction_smiles(reaction, reaction_format)
    elif "~" in any_smiles:
        # we have a multicomponent SMILES
        return apply_to_multicomponent_smiles(any_smiles, fn=fn, fragment_bond="~")
    else:
        # we have a single-component SMILES
        return fn(any_smiles)


def canonicalize_any(any_smiles: str, check_valence: bool = True) -> str:
    """
    Canonicalize any SMILES string (molecule SMILES, multicomponent SMILES, reaction SMILES).

    In the case of reaction SMILES, the format is kept.

    Args:
        any_smiles: any kind of SMILES string.
        check_valence: if False, will not do any valence check.

    Raises:
        Exception: different kinds of exception may be raised during parsing.
        InvalidSmiles: for canonicalization errors.

    Returns:
        the canonical (molecule, multicomponent, or reaction) SMILES string.
    """
    fn = partial(canonicalize_smiles, check_valence=check_valence)
    return apply_to_any_smiles(any_smiles, fn)


def sort_any(any_smiles: str) -> str:
    """
    Sort any SMILES string (molecule SMILES, multicomponent SMILES, reaction SMILES).

    For single-component SMILES, the fragments will be reorderd.
    In the case of reaction SMILES, the format is kept.

    Args:
        any_smiles: any kind of SMILES string.

    Raises:
        Exception: different kinds of exception may be raised during parsing.

    Returns:
        the sorted SMILES string.
    """
    if ">" in any_smiles:
        # we have a reaction SMILES
        reaction_format = determine_format(any_smiles)
        reaction = parse_reaction_smiles(any_smiles, reaction_format)
        reaction = sort_compounds(reaction)
        return to_reaction_smiles(reaction, reaction_format)
    else:
        # we call the same function for single- and multi-component SMILES
        return sort_multicomponent_smiles(any_smiles)
