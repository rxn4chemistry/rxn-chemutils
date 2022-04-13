from rxn_chemutils.utils import remove_atom_mapping


def test_remove_atom_mapping():
    # Reaction with mapping
    with_mapping = '[*:1]CCO.[CH3:1][C:2](=[O:3])[OH:4]>[H+]>CC[O:4][C:2](=[O:3])[CH3:1].O'
    without_mapping = '[*:1]CCO.[CH3][C](=[O])[OH]>[H+]>CC[O][C](=[O])[CH3].O'
    assert remove_atom_mapping(with_mapping) == without_mapping

    # Reaction with mapping -> do nothing
    rxn_without_mapping = '[*:1]CCO.CC(=O)O>[H+]>CCOC(C)=O.O'
    assert remove_atom_mapping(rxn_without_mapping) == rxn_without_mapping

    # Reaction with fragment info -> fragment info is not removed
    fragment_with_mapping = (
        '[*:1][Cl-].[Cl-].[Cl-].C(Cl)Cl>[CH3:9][CH:8]([CH3:10])[c:7]1[cH:11][cH:12][c:13]'
        '([cH:14][cH:15]1)[C:5](=[O:6])[CH2:4][CH2:3][CH2:2][Cl:1] |f:2.3.4.5|'
    )
    fragment_without_mapping = (
        '[*:1][Cl-].[Cl-].[Cl-].C(Cl)Cl>[CH3][CH]([CH3])[c]1[cH][cH][c]'
        '([cH][cH]1)[C](=[O])[CH2][CH2][CH2][Cl] |f:2.3.4.5|'
    )
    assert remove_atom_mapping(fragment_with_mapping) == fragment_without_mapping
