from typing import List

from rxn_chemutils.reaction_equation import (
    ReactionEquation, merge_reactants_and_agents, sort_compounds, canonicalize_compounds,
    remove_duplicate_compounds, has_repeated_molecules, rxn_standardization
)
from rxn_chemutils.conversion import canonicalize_smiles


def test_merge_reactants_and_agents():
    reactants = ['COCO', 'CCOC', 'OCC']
    agents = ['O', 'C', 'N']
    products = ['NCOC', 'CCNO', 'NCCC']

    reaction = ReactionEquation(reactants, agents, products)

    merged_reaction = merge_reactants_and_agents(reaction)

    expected = ReactionEquation(reactants + agents, [], products)
    assert merged_reaction == expected


def test_sort_compounds():
    reactants = ['COCO', 'CCOC', 'OCC']
    agents = ['O', 'C', 'N']
    products = ['NCOC', 'CCNO', 'NCCC']

    reaction = ReactionEquation(reactants, agents, products)

    sorted_reaction = sort_compounds(reaction)

    expected = ReactionEquation(sorted(reactants), sorted(agents), sorted(products))
    assert sorted_reaction == expected


def test_canonicalize_compounds():
    reactants = ['COCO', 'CCOC', 'OCC']
    agents = ['O', 'C', 'N']
    products = ['NCOC', 'CCNO', 'NCCC']

    reaction = ReactionEquation(reactants, agents, products)

    canonicalized_reaction = canonicalize_compounds(reaction)

    expected = ReactionEquation(
        [canonicalize_smiles(s) for s in reactants], [canonicalize_smiles(s) for s in agents],
        [canonicalize_smiles(s) for s in products]
    )
    assert canonicalized_reaction == expected


def test_remove_duplicate_compounds():
    reaction = ReactionEquation(
        ['A', 'B', 'C', 'A', 'B'], ['D', 'D', 'E'], ['A', 'B', 'G', 'H', 'H']
    )

    reaction_without_duplicates = remove_duplicate_compounds(reaction)
    # The compounds are in the same order, and duplicates in different categories are not removed
    expected = ReactionEquation(['A', 'B', 'C'], ['D', 'E'], ['A', 'B', 'G', 'H'])

    assert reaction_without_duplicates == expected


def test_equation_to_string():
    reactants = ['COCO', '[Na+].[OH-]', 'OCC']
    agents = ['O', 'C']
    products = ['NCOC']
    reaction = ReactionEquation(reactants, agents, products)

    reaction_string = reaction.to_string()
    expected = 'COCO.[Na+].[OH-].OCC>O.C>NCOC'

    assert reaction_string == expected


def test_equation_to_string_with_fragment_bond():
    reactants = ['COCO', '[Na+].[OH-]', 'OCC']
    agents = ['O', 'C']
    products = ['NCOC']
    reaction = ReactionEquation(reactants, agents, products)

    reaction_string = reaction.to_string(fragment_bond='~')
    expected = 'COCO.[Na+]~[OH-].OCC>O.C>NCOC'

    assert reaction_string == expected


def test_equation_from_string():
    reaction_string = 'COCO.[Na+].[OH-].OCC>O.C>NCOC'

    reaction = ReactionEquation.from_string(reaction_string)

    expected_reactants = ['COCO', '[Na+]', '[OH-]', 'OCC']
    expected_agents = ['O', 'C']
    expected_products = ['NCOC']
    expected_reaction = ReactionEquation(expected_reactants, expected_agents, expected_products)

    assert reaction == expected_reaction


def test_equation_from_string_with_fragment_bond():
    reaction_string = 'COCO.[Na+]~[OH-].OCC>O.C>NCOC'

    reaction = ReactionEquation.from_string(reaction_string, fragment_bond='~')

    expected_reactants = ['COCO', '[Na+].[OH-]', 'OCC']
    expected_agents = ['O', 'C']
    expected_products = ['NCOC']
    expected_reaction = ReactionEquation(expected_reactants, expected_agents, expected_products)

    assert reaction == expected_reaction


def test_equation_from_string_with_no_agent():
    reaction_string = 'COCO.[Na+].[OH-].OCC>>NCOC'

    reaction = ReactionEquation.from_string(reaction_string)

    expected_reactants = ['COCO', '[Na+]', '[OH-]', 'OCC']
    expected_agents: List[str] = []
    expected_products = ['NCOC']
    expected_reaction = ReactionEquation(expected_reactants, expected_agents, expected_products)

    assert reaction == expected_reaction


def test_has_repeated_molecules():
    assert not has_repeated_molecules(ReactionEquation(['A', 'B'], ['C'], ['D', 'E']))

    assert has_repeated_molecules(ReactionEquation(['A', 'A'], ['C'], ['D', 'E']))
    assert has_repeated_molecules(ReactionEquation(['A', 'B'], ['C', 'C'], ['D', 'E']))
    assert has_repeated_molecules(ReactionEquation(['A', 'B'], ['C'], ['D', 'E', 'E']))
    assert has_repeated_molecules(ReactionEquation(['A', 'B'], ['A'], ['D', 'E']))
    assert has_repeated_molecules(ReactionEquation(['A', 'B'], ['C'], ['B', 'E']))


def test_rxn_standardization():
    initial_reaction_equation = ReactionEquation(
        reactants=['OCC', 'O', 'C'],
        agents=['CCO'],
        products=['C'],
    )

    # NB:
    # 1. CCO is merged into the reactants
    # 2. OCC is canonicalized to CCO
    # 3. CCO is now present two times and removed
    # 4. sorted alphabetically
    expected = ReactionEquation(
        reactants=['C', 'CCO', 'O'],
        agents=[],
        products=['C'],
    )

    assert rxn_standardization(initial_reaction_equation) == expected
