from collections import Counter

from rxn_chemutils.miscellaneous import (equivalent_smiles, atom_type_counter, is_valid_smiles)


def test_equivalent_smiles():
    assert equivalent_smiles('CCO', 'C(C)O')
    assert not equivalent_smiles('CCO', 'O(C)C')

    assert equivalent_smiles('CCO', 'C(C)O', 'OCC')
    assert not equivalent_smiles('CCO', 'C(C)O', 'O(C)C')

    # Should be false for invalid SMILES
    assert not equivalent_smiles('CCO', '', 'O(C)C')


def test_is_valid_smiles():
    assert is_valid_smiles('CCO')
    assert is_valid_smiles('C')
    assert is_valid_smiles('NOOOOO')
    assert is_valid_smiles('CC(CCC)')

    assert not is_valid_smiles('YES')
    assert not is_valid_smiles('Noo')
    assert not is_valid_smiles('CC(CC')


def test_atom_type_counter():
    assert atom_type_counter('CCO') == Counter({'C': 2, 'H': 6, 'O': 1})
    assert atom_type_counter('CC=O') == Counter({'C': 2, 'H': 4, 'O': 1})
    assert atom_type_counter('[Na+].[Cl-]') == Counter({'Na': 1, 'Cl': 1})
