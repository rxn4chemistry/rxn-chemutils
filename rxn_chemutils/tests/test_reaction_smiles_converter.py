from rxn_chemutils.reaction_equation import ReactionEquation
from rxn_chemutils.reaction_smiles_converter import ReactionSmilesConverter


def test_from_reaction_smiles():
    reaction_smiles = '[CH3:7][C:4]1=[CH:3][C:2](=[C:1]([CH:6]=[CH:5]1)Br)[F:8]' \
                      '>[Li]CCCC.C1CCOC1' \
                      '>[CH3:7][C:4]1=[CH:3][C:2](=[C:1]([CH:6]=[CH:5]1)C(=O)O)[F:8]'

    reaction = ReactionSmilesConverter.from_reaction_smiles(reaction_smiles)

    assert reaction.reactants == ['CC1=CC(F)=C(Br)C=C1']
    assert reaction.agents == ['[Li]CCCC', 'C1CCOC1']
    assert reaction.products == ['CC1=CC(F)=C(C(=O)O)C=C1']


def test_from_reaction_smiles_with_fragments():
    reaction_smiles = 'BrC1=C(C(=O)O)C=C(C=C1)COC.C1(O)=CC(O)=CC=C1.[OH-].[Na+]' \
                      '>S(=O)(=O)([O-])[O-].[Cu+2].O' \
                      '>OC1=CC=C2C3=C(C(OC2=C1)=O)C=C(C=C3)COC' \
                      ' |f:2.3,4.5|'

    reaction = ReactionSmilesConverter.from_reaction_smiles(reaction_smiles)

    assert reaction.reactants == ['BrC1=C(C(=O)O)C=C(COC)C=C1', 'C1(O)=CC(O)=CC=C1', '[Na+].[OH-]']
    assert reaction.agents == ['O', 'O=S(=O)([O-])[O-].[Cu+2]']
    assert reaction.products == ['OC1=CC=C2C3=C(C(=O)OC2=C1)C=C(COC)C=C3']


def test_from_reaction_smiles_with_NaH():
    """
    In the original implementation using RDKit's ChemicalReaction, some molecules were parsed incorrectly, and this
    test would not have passed
    """
    reaction_smiles = 'C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3.N(=NC(=O)OCC)C(=O)OCC.C(CO)(=O)OCC.O=C1C(CCCC1)C#N.[H-].[Na+]' \
                      '>C1CCOC1.C1CCOC1.C1CCOC1.C1CCOC1' \
                      '>C(C)OC(=O)C=1OC2=C(C1N)CCCC2' \
                      ' |f:4.5|'

    reaction = ReactionSmilesConverter.from_reaction_smiles(reaction_smiles)

    assert reaction.reactants == [
        'C1=CC=C(P(C2=CC=CC=C2)C2=CC=CC=C2)C=C1', 'N(=NC(=O)OCC)C(=O)OCC', 'C(CO)(=O)OCC',
        'O=C1C(C#N)CCCC1', '[H-].[Na+]'
    ]
    assert reaction.agents == ['C1CCOC1', 'C1CCOC1', 'C1CCOC1', 'C1CCOC1']
    assert reaction.products == ['C(C)OC(=O)C1=C(N)C2=C(O1)CCCC2']


def test_to_reaction_smiles():
    reactants = ['CC1=CC(F)=C(Br)C=C1']
    agents = ['[Li]CCCC', 'C1CCOC1']
    products = ['CC1=CC(F)=C(C(=O)O)C=C1']

    reaction = ReactionEquation(reactants, agents, products)

    reaction_smiles = ReactionSmilesConverter.to_reaction_smiles(reaction)
    expected = 'CC1=CC(F)=C(Br)C=C1' \
               '>[Li]CCCC.C1CCOC1' \
               '>CC1=CC(F)=C(C(=O)O)C=C1'

    assert reaction_smiles == expected


def test_to_reaction_smiles_with_fragments():
    reactants = ['BrC1=C(C(=O)O)C=C(COC)C=C1', 'C1(O)=CC(O)=CC=C1', '[Na+].[OH-]']
    agents = ['O', 'O=S(=O)([O-])[O-].[Cu+2]']
    products = ['OC1=CC=C2C3=C(C(=O)OC2=C1)C=C(COC)C=C3']

    reaction = ReactionEquation(reactants, agents, products)

    reaction_smiles = ReactionSmilesConverter.to_reaction_smiles(reaction)
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

    reaction_smiles = ReactionSmilesConverter.to_reaction_smiles(reaction)
    expected = 'A.B.C.D.E' \
               '>F.G.H' \
               '>I.J.K.L.M.N' \
               ' |f:1.2,5.6,11.12.13|'

    assert reaction_smiles == expected
