import operator
import re
from functools import reduce
from typing import List, Optional, Tuple, Union

from rdkit import Chem, RDLogger
from rdkit.Chem import AssignStereochemistry, RemoveHs, SanitizeFlags, SanitizeMol
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import (
    MolFromMolBlock,
    MolFromSmiles,
    MolToMolBlock,
    MolToSmiles,
)

from .exceptions import InvalidMdl, InvalidSmiles, SanitizationError

RDLogger.logger().setLevel(RDLogger.CRITICAL)


def smiles_to_mol(
    smiles: str, sanitize: bool = True, find_radicals: bool = True
) -> Mol:
    """
    Convert a SMILES string to an RDKit Mol.

    Mainly a wrapper around MolFromSmiles that raises InvalidSmiles when necessary.

    Raises:
        InvalidSmiles for empty or invalid SMILES strings.

    Args:
        smiles: SMILES string to convert.
        sanitize: whether to sanitize the molecules or not. Note: sanitization here
            corresponds to doing a sanitization with SANITIZE_ALL.
        find_radicals: usually, it will be very practical to have this set to
            True, otherwise the imported Mol instances will potentially "add"
            hydrogen atoms to radical atoms. However, in some cases it may be
            useful to de-activate it because it may cause problems on aromatic
            molecules. Irrelevant if sanitize=True.

    Returns:
        Mol instance.
    """
    mol = MolFromSmiles(smiles, sanitize=sanitize)
    if not smiles or mol is None:
        raise InvalidSmiles(smiles)

    # MolFromSmiles with sanitize=False ignores all the radicals and converts
    # them back to normal atoms. To avoid this, we need to call the sanitization
    # function, either with no sanitization, or with radical finding.
    if not sanitize:
        if find_radicals:
            sanitizations = [Chem.SANITIZE_FINDRADICALS]
        else:
            sanitizations = [Chem.SANITIZE_NONE]

        sanitize_mol(mol, include_sanitizations=sanitizations)

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

    Mainly a wrapper around MolToMolBlock.
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
        raise ValueError(
            "Cannot specify both include_sanitizations and exclude_sanitizations."
        )

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
    mol = smiles_to_mol(smiles, sanitize=False, find_radicals=False)

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


def smiles_to_inchi(smiles: str, extended_tautomer_check: bool = False) -> str:
    """
    Get the InChI string for a given SMILES.

    Args:
        smiles: the SMILES to convert to InChI.
        extended_tautomer_check: include the options for additional tautomer standardization.

    Raises:
        InvalidSmiles for conversion errors or invalid SMILES.

    Returns:
        InChI string.
    """

    mol = smiles_to_mol(smiles, sanitize=False)

    # Due to a bug (?) in RDKit, it is necessary to reassign the stereochemistry
    # before conversion to InChi: https://github.com/rdkit/rdkit/issues/2361.
    AssignStereochemistry(mol, cleanIt=True, force=True)

    try:
        return mol_to_inchi(mol, extended_tautomer_check=extended_tautomer_check)
    except Exception:
        raise InvalidSmiles(smiles)


def mol_to_inchi(mol: Mol, extended_tautomer_check: bool = False) -> str:
    """
    Convert an RDKit Mol to an InChI.

    Args:
        mol: the RDKit Mol to convert to InChI.
        extended_tautomer_check: include the options for additional tautomer
            standardization.
    """

    options = ""
    if extended_tautomer_check:
        # https://pubs.acs.org/doi/10.1021/acs.jcim.9b01080
        options = "-KET -15T"

    return Chem.MolToInchi(mol, options=options)


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

    m = re.search(r"^(\S+) ?(.*)$", reaction_smiles)
    assert m is not None
    return m.group(1), m.group(2)
