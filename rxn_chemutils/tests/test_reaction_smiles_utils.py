from rxn_chemutils.reaction_smiles_utils import determine_fragment_groups, merge_molecules_from_fragment_groups


def test_fragment_groups():
    assert determine_fragment_groups('|f:2.3,4.5|') == [[2, 3], [4, 5]]
    assert determine_fragment_groups('|f:2.3,4.5.8|') == [[2, 3], [4, 5, 8]]
    assert determine_fragment_groups('|f:2.3,4|') == [[2, 3], [4]]
    assert determine_fragment_groups('|f:12.3,4.105|') == [[12, 3], [4, 105]]
    assert determine_fragment_groups('') == []


def test_merge_molecules_from_fragment_groups():
    # We consider, as an example, the following reaction
    # C.CC.CCC.CCCC > N.NN.NNN.NNNN.NNNNN > O.OO.OOO.OOOO
    # with fragmentation info |f:0.2.3,4.6,5.8|

    groups = [[0, 2, 3], [4, 6], [5, 8]]

    assert merge_molecules_from_fragment_groups(['C', 'CC', 'CCC', 'CCCC'], groups,
                                                0) == ['CC', 'C.CCC.CCCC']
    assert merge_molecules_from_fragment_groups(['N', 'NN', 'NNN', 'NNNN', 'NNNNN'], groups,
                                                4) == ['NNNN', 'N.NNN', 'NN.NNNNN']
    assert merge_molecules_from_fragment_groups(['O', 'OO', 'OOO', 'OOOO'], groups,
                                                9) == ['O', 'OO', 'OOO', 'OOOO']
