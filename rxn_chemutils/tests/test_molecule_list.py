from rxn_chemutils.molecule_list import MoleculeList


def test_molecule_list_behavior():
    molecules = ['A', 'B']

    # Equivalent to a list
    assert MoleculeList(molecules) == molecules

    # Random access
    assert MoleculeList(molecules)[0] == 'A'
    assert MoleculeList(molecules)[1] == 'B'

    # other checks
    assert list(reversed(MoleculeList(['A', 'B']))) == ['B', 'A']
    assert list(sorted(MoleculeList(['D', 'B', 'C', 'A']))) == ['A', 'B', 'C', 'D']


def test_string_conversion():
    assert MoleculeList.from_string('') == []
    assert MoleculeList.from_string('A.B.C') == ['A', 'B', 'C']
    assert MoleculeList.from_string('A~B.C') == ['A~B', 'C']
    assert MoleculeList.from_string('A~B.C', fragment_bond='~') == ['A.B', 'C']
    assert MoleculeList.from_string('A#B.C', fragment_bond='#') == ['A.B', 'C']

    assert MoleculeList([]).to_string() == ''
    assert MoleculeList(['A', 'B', 'C']).to_string() == 'A.B.C'
    assert MoleculeList(['A.B', 'C']).to_string() == 'A.B.C'
    assert MoleculeList(['A.B', 'C']).to_string(fragment_bond='~') == 'A~B.C'
    assert MoleculeList(['A.B', 'C']).to_string(fragment_bond='#') == 'A#B.C'


def test_to_list():
    molecule_list = MoleculeList(['A', 'B'])

    assert molecule_list == ['A', 'B']
    assert molecule_list.to_list() == ['A', 'B']

    assert type(molecule_list) == MoleculeList
    assert type(molecule_list.to_list()) == list
