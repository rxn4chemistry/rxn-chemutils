"""
Utilities related to "multi-component SMILES", i.e. strings containing multiple compounds
in SMILES notation, which may include fragment bonds.
"""
from functools import partial
from typing import Callable, Iterable, List, Optional

from rxn.utilities.containers import remove_duplicates

from .conversion import canonicalize_smiles


def multicomponent_smiles_to_list(
    multicomponent_smiles: str, fragment_bond: Optional[str] = None
) -> List[str]:
    """
    Convert a string of molecules into a list of molecules (taking fragment bonds into account).

    Args:
        multicomponent_smiles: multicomponent SMILES string to convert to a list.
        fragment_bond: fragment bond.

    Returns:
        The list of molecule SMILES comprised in the multi-component SMILES string.
    """
    molecules = multicomponent_smiles.split(".")
    molecules = [molecule for molecule in molecules if molecule != ""]

    # replace fragment bonds if necessary
    if fragment_bond is not None:
        molecules = [molecule.replace(fragment_bond, ".") for molecule in molecules]
    return molecules


def list_to_multicomponent_smiles(
    molecules: Iterable[str], fragment_bond: Optional[str] = None
) -> str:
    """
    Convert a list of molecules into a string representation (taking fragment
    bonds into account).

    Args:
        molecules: molecule SMILES strings to merge into a multi-component SMILES string.
        fragment_bond: fragment bond.

    Returns:
        A multi-component SMILES string.
    """
    # replace fragment bonds if necessary
    if fragment_bond is not None:
        molecules = [molecule.replace(".", fragment_bond) for molecule in molecules]

    return ".".join(molecules)


def apply_to_multicomponent_smiles(
    multicomponent_smiles: str,
    fn: Callable[[str], str],
    fragment_bond: Optional[str] = None,
) -> str:
    """
    Apply a function to the individual compounds in a multi-component SMILES string.

    Args:
        multicomponent_smiles: multi-component SMILES string to apply the function to.
        fn: function to apply on the distinct molecule SMILES.
        fragment_bond: fragment bond to use when parsing.

    Returns:
        New multi-component SMILES string after application of the function to the molecules.
    """
    molecules = multicomponent_smiles_to_list(
        multicomponent_smiles, fragment_bond=fragment_bond
    )
    molecules = [fn(molecule) for molecule in molecules]
    return list_to_multicomponent_smiles(molecules, fragment_bond=fragment_bond)


def canonicalize_multicomponent_smiles(
    multicomponent_smiles: str,
    fragment_bond: Optional[str] = None,
    check_valence: bool = True,
) -> str:
    """
    Canonicalize the molecules of a multi-component SMILES string.
    """
    canonicalize_fn = partial(canonicalize_smiles, check_valence=check_valence)
    return apply_to_multicomponent_smiles(
        multicomponent_smiles, canonicalize_fn, fragment_bond=fragment_bond
    )


def sort_multicomponent_smiles(multicomponent_smiles: str) -> str:
    """
    Sort the molecule SMILES in a multi-component SMILES string alphabetically.

    Note: no fragment bond is needed here, as it would not have any effect at all.
    """
    return list_to_multicomponent_smiles(
        sorted(multicomponent_smiles_to_list(multicomponent_smiles))
    )


def remove_duplicates_in_multicomponent_smiles(multicomponent_smiles: str) -> str:
    """
    Remove duplicate molecule SMILES strings in a multi-component SMILES string.

    Note: no fragment bond is needed here, as it would not have any effect at all.
    """
    return list_to_multicomponent_smiles(
        remove_duplicates(multicomponent_smiles_to_list(multicomponent_smiles))
    )
