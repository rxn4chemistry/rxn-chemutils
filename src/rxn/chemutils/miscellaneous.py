import re
import typing
from collections import Counter
from typing import List

from rdkit.Chem import AddHs, Atom, Mol

from .conversion import canonicalize_smiles, smiles_to_mol
from .exceptions import InvalidSmiles
from .multicomponent_smiles import canonicalize_multicomponent_smiles
from .reaction_equation import canonicalize_compounds
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


def canonicalize_any(any_smiles: str, check_valence: bool = True) -> str:
    """
    Canonicalize any SMILES string (molecule SMILES, multicomponent SMILES, reaction SMILES).

    In the case of reaction SMILES, the format is kept.

    Args:
        any_smiles: any kind of SMILES string.
        check_valence: if False, will not do any valence check.

    Raises:
        InvalidSmiles: if the SMILES string is not valid.

    Returns:
        the canonical (molecule, multicomponent, or reaction) SMILES string.
    """
    if ">" in any_smiles:
        # we have a reaction SMILES
        reaction_format = determine_format(any_smiles)
        reaction = parse_reaction_smiles(any_smiles, reaction_format)
        reaction = canonicalize_compounds(reaction, check_valence=check_valence)
        return to_reaction_smiles(reaction, reaction_format)
    else:
        # This covers both single molecule SMILES and multicomponent SMILES
        return canonicalize_multicomponent_smiles(
            any_smiles, fragment_bond="~", check_valence=check_valence
        )
