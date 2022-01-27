# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

from collections import Counter

import pytest

from rxn_chemutils.exceptions import InvalidSmiles
from rxn_chemutils.miscellaneous import (
    equivalent_smiles, atom_type_counter, is_valid_smiles, remove_atom_mapping,
    remove_chiral_centers
)


def test_equivalent_smiles():
    assert equivalent_smiles('CCO', 'C(C)O')
    assert not equivalent_smiles('CCO', 'O(C)C')

    assert equivalent_smiles('CCO', 'C(C)O', 'OCC')
    assert not equivalent_smiles('CCO', 'C(C)O', 'O(C)C')
    assert not equivalent_smiles('CF=C', 'F(C)C')

    # Should be false for invalid SMILES
    assert not equivalent_smiles('CCO', '', 'O(C)C')

    # In case of invalid valence: not equivalent if check_valence flag is specified
    assert equivalent_smiles('CFC', 'F(C)C')
    assert not equivalent_smiles('CFC', 'F(C)C', check_valence=True)


def test_is_valid_smiles():
    assert is_valid_smiles('CCO')
    assert is_valid_smiles('C')
    assert is_valid_smiles('NOOOOO')
    assert is_valid_smiles('CC(CCC)')

    assert not is_valid_smiles('CFC')
    assert not is_valid_smiles('YES')
    assert not is_valid_smiles('Noo')
    assert not is_valid_smiles('CC(CC')

    assert is_valid_smiles('CFC', check_valence=False)
    assert not is_valid_smiles('YES', check_valence=False)
    assert not is_valid_smiles('CC(CC', check_valence=False)


def test_atom_type_counter():
    assert atom_type_counter('CCO') == Counter({'C': 2, 'H': 6, 'O': 1})
    assert atom_type_counter('CC=O') == Counter({'C': 2, 'H': 4, 'O': 1})
    assert atom_type_counter('[Na+].[Cl-]') == Counter({'Na': 1, 'Cl': 1})

    # with radicals
    assert atom_type_counter('[C]CO') == Counter({'C': 2, 'H': 3, 'O': 1})

    with pytest.raises(InvalidSmiles):
        _ = atom_type_counter('CC(')


def test_remove_atom_mapping():
    # Reaction with mapping
    with_mapping = 'CCO.[CH3:1][C:2](=[O:3])[OH:4]>[H+]>CC[O:4][C:2](=[O:3])[CH3:1].O'
    without_mapping = 'CCO.[CH3][C](=[O])[OH]>[H+]>CC[O][C](=[O])[CH3].O'
    assert remove_atom_mapping(with_mapping) == without_mapping

    # Reaction with mapping -> do nothing
    rxn_without_mapping = 'CCO.CC(=O)O>[H+]>CCOC(C)=O.O'
    assert remove_atom_mapping(rxn_without_mapping) == rxn_without_mapping

    # Reaction with fragment info -> fragment info is not removed
    fragment_with_mapping = (
        '[Cl-].[Cl-].[Cl-].C(Cl)Cl>[CH3:9][CH:8]([CH3:10])[c:7]1[cH:11][cH:12][c:13]'
        '([cH:14][cH:15]1)[C:5](=[O:6])[CH2:4][CH2:3][CH2:2][Cl:1] |f:2.3.4.5|'
    )
    fragment_without_mapping = (
        '[Cl-].[Cl-].[Cl-].C(Cl)Cl>[CH3][CH]([CH3])[c]1[cH][cH][c]'
        '([cH][cH]1)[C](=[O])[CH2][CH2][CH2][Cl] |f:2.3.4.5|'
    )
    assert remove_atom_mapping(fragment_with_mapping) == fragment_without_mapping


def test_remove_chiral_centers():
    input_expected_dict = {
        'O[C@](Br)(C)N': 'O[C](Br)(C)N',
        'C[C@@](Br)(O)N': 'C[C](Br)(O)N',
        'N[C@H](O)C': 'N[CH](O)C',
        'N1[C@H](Cl)[C@@H](Cl)C(Cl)CC1': 'N1[CH](Cl)[CH](Cl)C(Cl)CC1',
        'F/C=C/F': 'F/C=C/F',
        '[N@+]': '[N+]',
        '[N@@+]': '[N+]',
        '[Si@@]': '[Si]',
        '[Si@H]': '[SiH]'
    }

    assert all(
        remove_chiral_centers(smi) == expected for smi, expected in input_expected_dict.items()
    )
