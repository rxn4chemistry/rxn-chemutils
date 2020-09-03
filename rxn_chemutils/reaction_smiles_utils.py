import re
from typing import List

from .utils import canonicalize_smiles


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
