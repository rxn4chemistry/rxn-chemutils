# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

import operator
import re
from functools import reduce
from typing import Sequence, List, Tuple, Optional, Union

from rdkit import RDLogger, Chem
from rdkit.Chem import (
    MolFromInchi, MolToInchi, SanitizeMol, SanitizeFlags, AssignStereochemistry, RemoveHs
)
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolToSmiles, MolFromSmiles, MolFromMolBlock, MolToMolBlock

from rxn_chemutils.exceptions import InvalidSmiles, SanitizationError, InvalidMdl

RDLogger.logger().setLevel(RDLogger.CRITICAL)


def smiles_to_mol(smiles: str, sanitize: bool = True) -> Mol:
    """
    Convert a SMILES string to an RDKit Mol.

    Mainly a wrapper around MolFromSmiles that raises InvalidSmiles when necessary.

    Raises:
        InvalidSmiles for empty or invalid SMILES strings.

    Args:
        smiles: SMILES string to convert.
        sanitize: whether to sanitize the molecules or not. Note: sanitization here
            corresponds to doing a sanitization with SANITIZE_ALL.

    Returns:
        Mol instance.
    """
    mol = MolFromSmiles(smiles, sanitize=sanitize)
    if not smiles or mol is None:
        raise InvalidSmiles(smiles)

    # MolFromSmiles with sanitize=False ignores all the radicals and converts
    # them back to normal atoms. To avoid this, we need to sanitize the radicals
    # manually for the case sanitize=False
    if not sanitize:
        sanitize_mol(mol, include_sanitizations=[Chem.SANITIZE_FINDRADICALS])

    return mol


def mol_to_smiles(mol: Mol, canonical: bool = True) -> str:
    """
    Convert an RDKit Mol to a SMILES string.

    Mainly a wrapper around MolToSmiles.
    """
    return MolToSmiles(mol, canonical=canonical)


def mdl_to_mol(mdl: str, sanitize: bool = True) -> Mol:
    """
    Convert an MDL Mol string to an RDKit Mol.

    Mainly a wrapper around MolFromMolBlock that raises InvalidMdl when necessary.

    Raises:
        InvalidMdl for empty or invalid MDL mol strings.

    Args:
        mdl: MDL Mol string to convert.
        sanitize: whether to sanitize the molecules or not. Note: sanitization here
            corresponds to doing a sanitization with SANITIZE_ALL.

    Returns:
        Mol instance.
    """
    mol = MolFromMolBlock(mdl, sanitize=sanitize)
    if mol is None:
        raise InvalidMdl(mdl)

    return mol


def mol_to_mdl(mol: Mol) -> str:
    """
    Convert an RDKit Mol to an MDL Mol string.

    Mainly a wrapper around MolToSmiles.
    """
    return MolToMolBlock(mol)


def sanitize_mol(
    mol: Mol,
    *,
    include_sanitizations: Optional[List[Union[SanitizeFlags, int]]] = None,
    exclude_sanitizations: Optional[List[Union[SanitizeFlags, int]]] = None
) -> None:
    """
    Sanitize an RDKit Mol with the specification of sanitizations to include or
    to exclude.

    Note: the RDKit sanitization function does not remove unnecessary hydrogens. See
    the function remove_hydrogens to do that.

    Raises:
        SanitizationError for unsuccessful sanitizations

    Args:
        mol: molecule to sanitize
        include_sanitizations: sanitizations to do. Exclusive with exclude_sanitizations.
        exclude_sanitizations: sanitizations to exclude, all the other ones will
            be applied. Exclusive with include_sanitizations.
    """

    # if no details about which sanitization steps to do is given, do them all
    if include_sanitizations is None and exclude_sanitizations is None:
        include_sanitizations = [Chem.SANITIZE_ALL]

    if include_sanitizations is not None and exclude_sanitizations is None:
        sanitize_ops = reduce(operator.or_, include_sanitizations, 0)
    elif include_sanitizations is None and exclude_sanitizations is not None:
        sanitize_ops = reduce(operator.xor, exclude_sanitizations, Chem.SANITIZE_ALL)
    else:
        raise ValueError('Cannot specify both include_sanitizations and exclude_sanitizations.')

    try:
        SanitizeMol(mol, sanitizeOps=sanitize_ops)
    except Exception as e:
        raise SanitizationError(mol) from e


def remove_hydrogens(mol: Mol) -> Mol:
    """
    Remove unnecessary hydrogens in a molecule.

    NB: The sanitization that is otherwise automatically done by RDKit is disabled.

    Returns:
        a new Mol object without unnecessary hydrogens.
    """
    return RemoveHs(mol, sanitize=False)


def canonicalize_smiles(smiles: str, check_valence: bool = True) -> str:
    """
    Canonicalize a SMILES string for a molecule.

    Args:
        smiles: SMILES string to canonicalize.
        check_valence: if False, will not do any valence check.

    Raises:
        InvalidSmiles for problems in parsing SMILES or in the sanitization.

    Returns:
        Canonicalized SMILES string.
    """
    mol = smiles_to_mol(smiles, sanitize=False)

    # NB: Removal of the unnecessary hydrogen atoms is disabled with sanitize=False above,
    # but the RDKit sanitize function does not actually do this. It is therefore
    # necessary to call this separately.
    mol = remove_hydrogens(mol)

    # Sanitization as a separate step, to enable exclusion of valence check
    try:
        excluded_sanitizations = []
        if not check_valence:
            excluded_sanitizations.append(Chem.SANITIZE_PROPERTIES)
        sanitize_mol(mol, exclude_sanitizations=excluded_sanitizations)
    except SanitizationError as e:
        raise InvalidSmiles(smiles) from e

    return mol_to_smiles(mol)


def maybe_canonicalize(smiles: str, check_valence: bool = True) -> str:
    """
    Canonicalize a SMILES string, but returns the original SMILES string if it fails.
    """
    try:
        return canonicalize_smiles(smiles, check_valence=check_valence)
    except InvalidSmiles:
        return smiles


def canonicalize_smiles_with_fragment_bonds(smiles: str, fragment_bond='~') -> str:
    """
    Canonicalizes a SMILES string that contains fragment bonds
    """
    # Raise for empty SMILES
    if not smiles:
        raise InvalidSmiles(smiles)

    try:
        return '.'.join(
            sorted(
                [
                    canonicalize_smiles(_.replace(fragment_bond, '.')).replace('.', fragment_bond)
                    for _ in smiles.split('.')
                ]
            )
        )
    except Exception:
        raise InvalidSmiles(smiles)


def smiles_to_inchi(smiles: str) -> str:
    """
    Get the InChI string for a given SMILES.

    Raises:
        InvalidSmiles for conversion errors or invalid SMILES.

    Returns:
        InChI string.
    """

    mol = smiles_to_mol(smiles, sanitize=False)

    # Due to a bug (?) in RDKit, it is necessary to reassign the stereochemistry
    # before conversion to InChi: https://github.com/rdkit/rdkit/issues/2361.
    # Strangely, it is also necessary to do an empty sanitization (NONE) beforehand.
    sanitize_mol(mol, include_sanitizations=[Chem.SANITIZE_NONE])
    AssignStereochemistry(mol, cleanIt=True, force=True)

    try:
        return MolToInchi(mol)
    except Exception:
        raise InvalidSmiles(smiles)


def inchify_smiles(smiles: str) -> str:
    """
    Inchify a SMILES string for a molecule
    """
    # Raise for empty SMILES
    if not smiles:
        raise InvalidSmiles(smiles)

    try:
        return MolToSmiles(MolFromInchi(MolToInchi(MolFromSmiles(smiles))))
    except Exception:
        raise InvalidSmiles(smiles)


def inchify_smiles_with_fragment_bonds(smiles: str, fragment_bond='~') -> str:
    """
    Inchify a SMILES string that contains fragment bonds for a molecule
    """
    # Raise for empty SMILES
    if not smiles:
        raise InvalidSmiles(smiles)

    try:
        return '.'.join(
            sorted(
                [
                    inchify_smiles(_.replace(fragment_bond, '.')).replace('.', fragment_bond)
                    for _ in smiles.split('.')
                ]
            )
        )
    except Exception:
        raise InvalidSmiles(smiles)


def cleanup_smiles(smiles: str) -> str:
    """
    Cleanup a SMILES string, doing the bare minimum.

    This means that no canonicalization, no valence check, no kekulization, etc,
    will be done.
    See the unit tests for examples.

    A minimal sanitization (SANITIZE_FINDRADICALS) is necessary, otherwise
    "[C]" is converted to "C".

    Args:
        smiles: SMILES to clean up.

    Returns:
        A cleaned-up SMILES string.
    """
    mol = smiles_to_mol(smiles, sanitize=False)
    sanitize_mol(mol, include_sanitizations=[Chem.SANITIZE_FINDRADICALS])
    return mol_to_smiles(mol, canonical=False)


def mols_to_smiles(mols: Sequence[Mol], canonical: bool = True) -> List[str]:
    return [mol_to_smiles(mol, canonical) for mol in mols]


def convert_group_to_mols(smiles_group: str) -> List[Mol]:
    """
    Args:
        smiles_group: List of SMILES strings, separated by dots
    """
    # Handle empty strings:
    if not smiles_group:
        return []

    smiles_list = smiles_group.split('.')
    return [smiles_to_mol(smiles) for smiles in smiles_list]


def split_smiles_and_fragment_info(reaction_smiles: str) -> Tuple[str, str]:
    """
    The reaction SMILES from Pistachio sometimes contain fraction information
    at the end of the given string. This function splits both parts of the
    reaction SMILES.

    Args:
        reaction_smiles: (potentially extended) reaction SMILES.

    Returns:
        Tuple: ('pure' reaction SMILES, fragment information).
    """

    m = re.search(r'^(\S+) ?(.*)$', reaction_smiles)
    assert m is not None
    return m.group(1), m.group(2)
