from rxn_chemutils.extended_reaction_smiles import (
    parse_extended_reaction_smiles, to_extended_reaction_smiles, determine_fragment_groups,
    merge_molecules_from_fragment_groups
)
from rxn_chemutils.reaction_equation import ReactionEquation


def test_from_reaction_smiles():
    reaction_smiles = '[CH3:7][C:4]1=[CH:3][C:2](=[C:1]([CH:6]=[CH:5]1)Br)[F:8]' \
                      '>[Li]CCCC.C1CCOC1' \
                      '>[CH3:7][C:4]1=[CH:3][C:2](=[C:1]([CH:6]=[CH:5]1)C(=O)O)[F:8]'

    reaction = parse_extended_reaction_smiles(reaction_smiles)

    assert reaction.reactants == ['CC1=CC(F)=C(Br)C=C1']
    assert reaction.agents == ['[Li]CCCC', 'C1CCOC1']
    assert reaction.products == ['CC1=CC(F)=C(C(=O)O)C=C1']


def test_from_reaction_smiles_with_fragments():
    reaction_smiles = 'BrC1=C(C(=O)O)C=C(C=C1)COC.C1(O)=CC(O)=CC=C1.[OH-].[Na+]' \
                      '>S(=O)(=O)([O-])[O-].[Cu+2].O' \
                      '>OC1=CC=C2C3=C(C(OC2=C1)=O)C=C(C=C3)COC' \
                      ' |f:2.3,4.5|'

    reaction = parse_extended_reaction_smiles(reaction_smiles)

    assert reaction.reactants == ['BrC1=C(C(=O)O)C=C(COC)C=C1', 'C1(O)=CC(O)=CC=C1', '[OH-].[Na+]']
    assert reaction.agents == ['O', 'S(=O)(=O)([O-])[O-].[Cu+2]']
    assert reaction.products == ['OC1=CC=C2C3=C(C(=O)OC2=C1)C=C(COC)C=C3']


def test_from_reaction_smiles_does_not_sanitize():
    # In earlier code, molecules were canonicalized if they were part of a fragment
    reaction_smiles = 'C1=CC=CC=C1.[N+](=O)(O)[O-]>>C1=CC=CC=C1N(=O)=O |f:0.1|'

    reaction = parse_extended_reaction_smiles(reaction_smiles)

    assert reaction.reactants == ['C1=CC=CC=C1.[N+](=O)(O)[O-]']
    assert reaction.agents == []
    assert reaction.products == ['C1=CC=CC=C1N(=O)=O']


def test_from_reaction_smiles_with_NaH():
    """
    In the original implementation using RDKit's ChemicalReaction, some molecules were parsed incorrectly, and this
    test would not have passed
    """
    reaction_smiles = 'C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3.N(=NC(=O)OCC)C(=O)OCC.C(CO)(=O)OCC.O=C1C(CCCC1)C#N.[H-].[Na+]' \
                      '>C1CCOC1.C1CCOC1.C1CCOC1.C1CCOC1' \
                      '>C(C)OC(=O)C=1OC2=C(C1N)CCCC2' \
                      ' |f:4.5|'

    reaction = parse_extended_reaction_smiles(reaction_smiles)

    assert reaction.reactants == [
        'C1=CC=C(P(C2=CC=CC=C2)C2=CC=CC=C2)C=C1', 'N(=NC(=O)OCC)C(=O)OCC', 'C(CO)(=O)OCC',
        'O=C1C(C#N)CCCC1', '[H-].[Na+]'
    ]
    assert reaction.agents == ['C1CCOC1', 'C1CCOC1', 'C1CCOC1', 'C1CCOC1']
    assert reaction.products == ['C(C)OC(=O)C1=C(N)C2=C(O1)CCCC2']


def test_from_reaction_smiles_with_rdkit_problem():
    # A previous version relying more on RDKit returned, as a reactant, 'C.[Pd]'.
    reaction_smiles = '[C].[Pd]>CC>CCO |f:0.1|'
    reaction = parse_extended_reaction_smiles(reaction_smiles)
    assert reaction.reactants == ['[C].[Pd]']
    assert reaction.agents == ['CC']
    assert reaction.products == ['CCO']

    # Related example from Pistachio
    reaction_smiles = (
        '[CH3:24][C:23]([CH3:25])([CH3:26])[O:22][C:20](=[O:21])[N:10]1[C@@H:11]([CH2:12][CH2:13]'
        '[C@@H:9]1[C:4]2=[CH:3][C:2](=[C:7]([CH:6]=[CH:5]2)[F:8])[F:1])/[CH:14]=[CH:15]/[C:16]'
        '(=[O:17])[O:18][CH3:19]>[HH].[C].CCOC(=O)C.[Pd]>[CH3:24][C:23]([CH3:25])([CH3:26])[O:22]'
        '[C:20](=[O:21])[N:10]1[C@@H:11]([CH2:12][CH2:13][C@@H:9]1[C:4]2=[CH:3][C:2](=[C:7]([CH:6]'
        '=[CH:5]2)[F:8])[F:1])[CH2:14][CH2:15][C:16](=[O:17])[O:18][CH3:19] |f:2.4|'
    )
    reaction = parse_extended_reaction_smiles(reaction_smiles)
    assert reaction.reactants == ['CC(C)(C)OC(=O)N1[C@H](/C=C/C(=O)OC)CC[C@@H]1C1=CC(F)=C(F)C=C1']
    assert reaction.agents == ['[HH]', 'CCOC(=O)C', '[C].[Pd]']
    assert reaction.products == ['CC(C)(C)OC(=O)N1[C@H](CCC(=O)OC)CC[C@@H]1C1=CC(F)=C(F)C=C1']


def test_to_reaction_smiles():
    reactants = ['CC1=CC(F)=C(Br)C=C1']
    agents = ['[Li]CCCC', 'C1CCOC1']
    products = ['CC1=CC(F)=C(C(=O)O)C=C1']

    reaction = ReactionEquation(reactants, agents, products)

    reaction_smiles = to_extended_reaction_smiles(reaction)
    expected = 'CC1=CC(F)=C(Br)C=C1' \
               '>[Li]CCCC.C1CCOC1' \
               '>CC1=CC(F)=C(C(=O)O)C=C1'

    assert reaction_smiles == expected


def test_to_reaction_smiles_with_fragments():
    reactants = ['BrC1=C(C(=O)O)C=C(COC)C=C1', 'C1(O)=CC(O)=CC=C1', '[Na+].[OH-]']
    agents = ['O', 'O=S(=O)([O-])[O-].[Cu+2]']
    products = ['OC1=CC=C2C3=C(C(=O)OC2=C1)C=C(COC)C=C3']

    reaction = ReactionEquation(reactants, agents, products)

    reaction_smiles = to_extended_reaction_smiles(reaction)
    expected = 'BrC1=C(C(=O)O)C=C(COC)C=C1.C1(O)=CC(O)=CC=C1.[Na+].[OH-]' \
               '>O.O=S(=O)([O-])[O-].[Cu+2]' \
               '>OC1=CC=C2C3=C(C(=O)OC2=C1)C=C(COC)C=C3' \
               ' |f:2.3,5.6|'

    assert reaction_smiles == expected


def test_to_reaction_smiles_with_dummy_fragments():
    reactants = ['A', 'B.C', 'D', 'E']
    agents = ['F.G', 'H']
    products = ['I', 'J', 'K', 'L.M.N']

    reaction = ReactionEquation(reactants, agents, products)

    reaction_smiles = to_extended_reaction_smiles(reaction)
    expected = 'A.B.C.D.E' \
               '>F.G.H' \
               '>I.J.K.L.M.N' \
               ' |f:1.2,5.6,11.12.13|'

    assert reaction_smiles == expected


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
