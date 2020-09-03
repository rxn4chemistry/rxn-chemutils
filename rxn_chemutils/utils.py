import re
import typing
from collections import Counter
from typing import List, Sequence, Tuple

from rdkit import RDLogger
from rdkit.Chem import MolFromSmiles, Mol, MolToSmiles, MolToInchi, MolFromInchi, Atom, AddHs

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


def remove_atom_mapping(mols: Sequence[Mol]) -> None:
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                atom.ClearProp('molAtomMapNumber')


def split_smiles_and_fragment_info(reaction_smiles) -> Tuple[str, str]:
    """
    The reaction SMILES from Pistachio sometimes contain fraction information at the end of the given string.
    This function splits both parts of the reaction SMILES.

    Returns: Tuple ('pure' reaction SMILES, fragment information)
    """

    m = re.search(r'^(\S+) ?(.*)$', reaction_smiles)
    assert m is not None
    return m.group(1), m.group(2)


def tokenize_smiles(smiles: str) -> List[str]:
    """
    Tokenize a SMILES molecule or reaction
    """
    pattern = r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smiles)]

    assert smiles == ''.join(tokens), f'SMI: {smiles} =!= JOIN: {"".join(tokens)}'

    return tokens


def smi_tokenizer(smi: str) -> str:
    """
    Tokenize a SMILES molecule or reaction, and join the tokens with spaces.
    """
    tokens = tokenize_smiles(smi)
    return ' '.join(tokens)


def smi_detokenizer(tokenized_smiles: str) -> str:
    """
    Detokenize a tokenized SMILES string (that contains spaces between the characters).

    Args:
        tokenized_smiles: tokenized SMILES, for instance 'C C ( C O ) = N >> C C ( C = O ) N'

    Returns:
        SMILES after detokenization, for instance 'CC(CO)=N>>CC(C=O)N'
    """
    return tokenized_smiles.replace(' ', '')


def atom_type_counter(smiles: str) -> typing.Counter[str]:
    """
    Returns a counter of atom types (as symbols).
    """

    mol: Mol = AddHs(MolFromSmiles(smiles))
    atoms: List[Atom] = mol.GetAtoms()
    return Counter(atom.GetSymbol() for atom in atoms)
