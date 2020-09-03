import re
from typing import Sequence, List, Tuple

from rdkit import RDLogger
from rdkit.Chem import MolFromInchi, MolToInchi
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolToSmiles, MolFromSmiles

from rxn_chemutils.rdkit_utils import remove_atom_mapping

RDLogger.logger().setLevel(RDLogger.CRITICAL)


class InvalidSmiles(ValueError):

    def __init__(self, smiles: str):
        super().__init__(f'"{smiles}" is not a valid SMILES string')


def canonicalize_smiles(smiles: str) -> str:
    """
    Canonicalizes a SMILES string for a molecule
    """
    # Raise for empty SMILES
    if not smiles:
        raise InvalidSmiles(smiles)

    try:
        return MolToSmiles(MolFromSmiles(smiles))
    except Exception:
        raise InvalidSmiles(smiles)


def maybe_canonicalize(smiles: str) -> str:
    """
    Canonicalize a SMILES string, but returns the original SMILES string if it fails.
    """
    try:
        return canonicalize_smiles(smiles)
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


def canonicalize_reaction_smiles(reaction_smiles: str) -> str:
    """
    Canonicalizes a SMILES string for a reaction.
    This will remove the atom mapping, but keep the fragment information.
    """
    reaction_smiles, fragment_info = split_smiles_and_fragment_info(reaction_smiles)

    raw_smiles_groups = reaction_smiles.split('>')
    assert len(raw_smiles_groups) == 3

    mol_groups = [convert_group_to_mols(smiles_group) for smiles_group in raw_smiles_groups]
    for group in mol_groups:
        remove_atom_mapping(group)

    smiles_groups = [mols_to_smiles(group) for group in mol_groups]

    joined_strings = ['.'.join(group) for group in smiles_groups]

    canonical_reaction_smiles = '>'.join(joined_strings)

    if fragment_info:
        canonical_reaction_smiles += f' {fragment_info}'

    return canonical_reaction_smiles


def smiles_to_mol(smiles: str, sanitize: bool = True) -> Mol:
    # Raise for empty SMILES
    if not smiles:
        raise InvalidSmiles(smiles)

    mol = MolFromSmiles(smiles, sanitize=sanitize)
    if mol is None:
        raise InvalidSmiles(smiles)

    return mol


def mol_to_smiles(mol: Mol, canonical: bool = True) -> str:
    return MolToSmiles(mol, canonical=canonical)


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
