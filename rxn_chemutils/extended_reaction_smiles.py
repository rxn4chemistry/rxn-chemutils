import re
from typing import List, Tuple

from rdkit.Chem import Mol

from .chemical_reaction import ChemicalReaction
from .conversion import mols_to_smiles, canonicalize_smiles, split_smiles_and_fragment_info
from .rdkit_utils import remove_atom_mapping
from .reaction_equation import ReactionEquation


class UnsupportedExtendedReactionSmiles(ValueError):

    def __init__(self, reaction_smiles: str):
        super().__init__(f'The syntax of "{reaction_smiles}" is not supported by RDKit.')


def parse_extended_reaction_smiles(extended_reaction_smiles: str) -> ReactionEquation:
    """
    Convert an extended reaction SMILES (with potential fragment information)
    to a ReactionEquation instance.

    Args:
        extended_reaction_smiles: extended reaction SMILES

    Returns:
        ReactionEquation instance
    """
    return _Importer.convert(extended_reaction_smiles)


def to_extended_reaction_smiles(reaction: ReactionEquation) -> str:
    """
    Convert a ReactionEquation instance to an extended reaction SMILES (with
    potential fragment information).

    Args:
        reaction: reaction equation to convert

    Returns:
        The extended reaction SMILES string.
    """
    return _Exporter.convert(reaction)


class _Importer:
    """
    Convert extended reaction SMILES to ReactionEquation instances.
    """

    @staticmethod
    def convert(reaction_smiles: str) -> ReactionEquation:
        pure_smiles, fragment_info = split_smiles_and_fragment_info(reaction_smiles)

        # Parse the reaction with RDKit
        try:
            rxn = ChemicalReaction(pure_smiles, sanitize=False)
        except Exception:
            raise UnsupportedExtendedReactionSmiles(reaction_smiles)

        try:
            return _Importer.process_reaction_participants(rxn, fragment_info)
        except Exception:
            raise UnsupportedExtendedReactionSmiles(reaction_smiles)

    @staticmethod
    def process_reaction_participants(
        rxn: ChemicalReaction, fragment_info: str
    ) -> ReactionEquation:
        raw_mol_groups = [rxn.reactants, rxn.agents, rxn.products]
        raw_smiles_groups = [_Importer.convert_to_smiles(group) for group in raw_mol_groups]

        fragment_groups = determine_fragment_groups(fragment_info)
        groups = _Importer.group_fragments(raw_smiles_groups, fragment_groups)

        return ReactionEquation(*groups)

    @staticmethod
    def convert_to_smiles(mols: List[Mol]) -> List[str]:
        remove_atom_mapping(mols)
        return mols_to_smiles(mols, canonical=False)

    @staticmethod
    def group_fragments(raw_smiles_groups: List[List[str]],
                        fragment_groups: List[List[int]]) -> List[List[str]]:
        """
        Merge the reaction fragments belonging together.
        """
        merged_groups = []

        offset = 0
        for raw_smiles_group in raw_smiles_groups:
            merged_smiles_group = merge_molecules_from_fragment_groups(
                raw_smiles_group, fragment_groups, offset
            )
            merged_groups.append(merged_smiles_group)
            offset += len(raw_smiles_group)

        return merged_groups


class _Exporter:
    """
    Convert ReactionEquation to extended reaction SMILES with fragment information.
    """

    @staticmethod
    def convert(reaction: ReactionEquation) -> str:
        offset = 0
        reactants, reactant_groups = _Exporter.fragment_group(reaction.reactants, offset)
        offset += len(reactants)
        agents, agent_groups = _Exporter.fragment_group(reaction.agents, offset)
        offset += len(agents)
        products, product_groups = _Exporter.fragment_group(reaction.products, offset)

        groups = reactant_groups + agent_groups + product_groups

        smiles_groups = (
            '.'.join(smiles for smiles in group) for group in (reactants, agents, products)
        )
        smiles_without_fragment_info = '>'.join(smiles_groups)

        if not groups:
            return smiles_without_fragment_info

        fragment_info = _Exporter.generate_fragment_info(groups)
        return f'{smiles_without_fragment_info} {fragment_info}'

    @staticmethod
    def fragment_group(compounds: List[str], offset: int) -> Tuple[List[str], List[List[int]]]:
        """
        Converts a group of molecules, some of which possibly are composed of several fragments,
        to the list of SMILES to put in the final reaction SMILES, and a list of groups of molecules
        belonging together.

        Example:
            ['O', '[Na+].[OH-]'] -> (['O', '[Na+]', '[OH-]'], [[1, 2]])

        Args:
            compounds: SMILES string for molecules of the same group (f.i. reactants)
            offset: index, in the final reaction SMILES, of the first molecule of that group
        """

        smiles_list: List[str] = []
        groups: List[List[int]] = []

        current_index = offset
        for c in compounds:
            molecules = c.split('.')
            smiles_list.extend(molecules)

            number_fragments = len(molecules)
            if number_fragments > 1:
                groups.append(list(range(current_index, current_index + number_fragments)))

            current_index += number_fragments

        return smiles_list, groups

    @staticmethod
    def generate_fragment_info(groups: List[List[int]]) -> str:
        if not groups:
            return ''

        group_strings = ['.'.join(str(number) for number in g) for g in groups]
        all_groups = ','.join(group_strings)
        return f'|f:{all_groups}|'


def determine_fragment_groups(fragment_info: str) -> List[List[int]]:
    """
    From the fragment info string (such as '|f:0.2,5.6|'), determine the groups of indices that must be grouped.

    Returns:
        List of groups (f.i. [[0,2], [5,6]])
    """

    m = re.findall(r'(\d+(?:\.\d+)*)', fragment_info)

    groups = []
    for match in m:
        indices = match.split('.')
        groups.append([int(i) for i in indices])
    return groups


def merge_molecules_from_fragment_groups(
    smiles_list: List[str], fragment_groups: List[List[int]], offset: int
) -> List[str]:
    """
    Combine molecules according to the fragment definition.

    Args:
        smiles_list: SMILES strings to potentially merge
        fragment_groups: indices for groups of molecules belonging together; some may be out of range
        offset: what index the first molecule in smiles_list corresponds to

    Returns:
        List of molecules, some of which were potentially merged. Non-merged molecules come first.
    """

    allowed_indices = set(range(len(smiles_list)))
    merged_indices = set()

    merged_molecules = []

    for group in fragment_groups:
        relative_indices = [i - offset for i in group]
        in_range = [e in allowed_indices for e in relative_indices]

        # either all indices are in the range, or none is
        all_in_range = all(in_range)
        none_in_range = not any(in_range)
        if not (all_in_range or none_in_range):
            raise ValueError()

        if all_in_range:
            merged_molecule = '.'.join(smiles_list[i] for i in relative_indices)
            merged_molecules.append(canonicalize_smiles(merged_molecule))
            merged_indices.update(relative_indices)

    remaining_molecule_indices = sorted(list(allowed_indices - merged_indices))
    unmerged_molecules = [smiles_list[i] for i in remaining_molecule_indices]

    return unmerged_molecules + merged_molecules
