from typing import List

import pytest

from rxn.chemutils.conversion import canonicalize_smiles
from rxn.chemutils.exceptions import InvalidReactionSmiles, InvalidSmiles
from rxn.chemutils.reaction_equation import (
    ReactionEquation,
    apply_to_compound_groups,
    apply_to_compounds,
    canonicalize_compounds,
    cleanup_compounds,
    has_repeated_molecules,
    merge_reactants_and_agents,
    remove_duplicate_compounds,
    remove_precursors_from_products,
    rxn_standardization,
    sort_compounds,
)


def test_merge_reactants_and_agents() -> None:
    reactants = ["COCO", "CCOC", "OCC"]
    agents = ["O", "C", "N"]
    products = ["NCOC", "CCNO", "NCCC"]

    reaction = ReactionEquation(reactants, agents, products)

    merged_reaction = merge_reactants_and_agents(reaction)

    expected = ReactionEquation(reactants + agents, [], products)
    assert merged_reaction == expected


def test_sort_compounds() -> None:
    reactants = ["COCO", "CCOC", "OCC"]
    agents = ["O", "C", "N"]
    products = ["NCOC", "CCNO", "NCCC"]

    reaction = ReactionEquation(reactants, agents, products)

    sorted_reaction = sort_compounds(reaction)

    expected = ReactionEquation(sorted(reactants), sorted(agents), sorted(products))
    assert sorted_reaction == expected


def test_canonicalize_compounds() -> None:
    reactants = ["COCO", "CCOC", "OCC"]
    agents = ["O", "C", "N"]
    products = ["NCOC", "CCNO", "NCCC"]

    reaction = ReactionEquation(reactants, agents, products)

    canonicalized_reaction = canonicalize_compounds(reaction)

    expected = ReactionEquation(
        [canonicalize_smiles(s) for s in reactants],
        [canonicalize_smiles(s) for s in agents],
        [canonicalize_smiles(s) for s in products],
    )
    assert canonicalized_reaction == expected


def test_canonicalize_compounds_with_invalid_valence() -> None:
    reaction = ReactionEquation(reactants=["F", "CC"], agents=[], products=["CFC"])

    with pytest.raises(InvalidSmiles):
        canonicalize_compounds(reaction)

    # does not raise if no valence check, returns original reaction
    assert canonicalize_compounds(reaction, check_valence=False) == reaction


def test_cleanup_compounds() -> None:
    reaction = ReactionEquation(["CC[CH2]C", "C-C"], [], ["[C].[Pd]"])

    clean_reaction = cleanup_compounds(reaction)
    expected = ReactionEquation(["CCCC", "CC"], [], ["[C].[Pd]"])

    assert clean_reaction == expected


def test_remove_duplicate_compounds() -> None:
    reaction = ReactionEquation(
        ["A", "B", "C", "A", "B"], ["D", "D", "E"], ["A", "B", "G", "H", "H"]
    )

    reaction_without_duplicates = remove_duplicate_compounds(reaction)
    # The compounds are in the same order, and duplicates in different categories are not removed
    expected = ReactionEquation(["A", "B", "C"], ["D", "E"], ["A", "B", "G", "H"])

    assert reaction_without_duplicates == expected


def test_equation_to_string() -> None:
    reactants = ["COCO", "[Na+].[OH-]", "OCC"]
    agents = ["O", "C"]
    products = ["NCOC"]
    reaction = ReactionEquation(reactants, agents, products)

    reaction_string = reaction.to_string()
    expected = "COCO.[Na+].[OH-].OCC>O.C>NCOC"

    assert reaction_string == expected


def test_equation_to_string_with_fragment_bond() -> None:
    reactants = ["COCO", "[Na+].[OH-]", "OCC"]
    agents = ["O", "C"]
    products = ["NCOC"]
    reaction = ReactionEquation(reactants, agents, products)

    reaction_string = reaction.to_string(fragment_bond="~")
    expected = "COCO.[Na+]~[OH-].OCC>O.C>NCOC"

    assert reaction_string == expected


def test_equation_from_string() -> None:
    reaction_string = "COCO.[Na+].[OH-].OCC>O.C>NCOC"

    reaction = ReactionEquation.from_string(reaction_string)

    expected_reactants = ["COCO", "[Na+]", "[OH-]", "OCC"]
    expected_agents = ["O", "C"]
    expected_products = ["NCOC"]
    expected_reaction = ReactionEquation(
        expected_reactants, expected_agents, expected_products
    )

    assert reaction == expected_reaction


def test_equation_from_string_ignores_additional_dots() -> None:
    reaction_smiles = "..A.B.>.>C.D."
    expected = ReactionEquation(["A", "B"], [], ["C", "D"])
    assert ReactionEquation.from_string(reaction_smiles) == expected


def test_equation_from_string_with_invalid_strings() -> None:
    invalid_reaction_smiles = ["A.B>C", "A>>>C", "A>>B>>C", "A.B.C"]
    for reaction_smiles in invalid_reaction_smiles:
        with pytest.raises(InvalidReactionSmiles):
            _ = ReactionEquation.from_string(reaction_smiles)


def test_equation_from_string_with_fragment_bond() -> None:
    reaction_string = "COCO.[Na+]~[OH-].OCC>O.C>NCOC"

    reaction = ReactionEquation.from_string(reaction_string, fragment_bond="~")

    expected_reactants = ["COCO", "[Na+].[OH-]", "OCC"]
    expected_agents = ["O", "C"]
    expected_products = ["NCOC"]
    expected_reaction = ReactionEquation(
        expected_reactants, expected_agents, expected_products
    )

    assert reaction == expected_reaction


def test_equation_from_string_with_no_agent() -> None:
    reaction_string = "COCO.[Na+].[OH-].OCC>>NCOC"

    reaction = ReactionEquation.from_string(reaction_string)

    expected_reactants = ["COCO", "[Na+]", "[OH-]", "OCC"]
    expected_agents: List[str] = []
    expected_products = ["NCOC"]
    expected_reaction = ReactionEquation(
        expected_reactants, expected_agents, expected_products
    )

    assert reaction == expected_reaction


def test_equation_from_string_with_dative_bond() -> None:
    reaction_string = "COC(=O)CCBr.O=C([O-]->[K+])>>COC(=O)CCOC(=O)C"

    reaction = ReactionEquation.from_string(reaction_string)

    expected_reactants = ["COC(=O)CCBr", "O=C([O-]->[K+])"]
    expected_agents: List[str] = []
    expected_products = ["COC(=O)CCOC(=O)C"]
    expected_reaction = ReactionEquation(
        expected_reactants, expected_agents, expected_products
    )

    assert reaction == expected_reaction


def test_iter_all_smiles() -> None:
    reaction_string = "COCO.[Na+].[OH-]>O>NCOC"
    reaction = ReactionEquation.from_string(reaction_string)
    assert list(reaction.iter_all_smiles()) == ["COCO", "[Na+]", "[OH-]", "O", "NCOC"]


def test_has_repeated_molecules() -> None:
    assert not has_repeated_molecules(ReactionEquation(["A", "B"], ["C"], ["D", "E"]))

    assert has_repeated_molecules(ReactionEquation(["A", "A"], ["C"], ["D", "E"]))
    assert has_repeated_molecules(ReactionEquation(["A", "B"], ["C", "C"], ["D", "E"]))
    assert has_repeated_molecules(ReactionEquation(["A", "B"], ["C"], ["D", "E", "E"]))
    assert has_repeated_molecules(ReactionEquation(["A", "B"], ["A"], ["D", "E"]))
    assert has_repeated_molecules(ReactionEquation(["A", "B"], ["C"], ["B", "E"]))


def test_apply_to_compounds() -> None:
    def dummy_fn(smiles: str) -> str:
        # add a 'C' to the SMILES
        return smiles + "C"

    assert apply_to_compounds(
        ReactionEquation(["C", "N", "O"], ["P"], ["CNO"]), dummy_fn
    ) == ReactionEquation(["CC", "NC", "OC"], ["PC"], ["CNOC"])


def test_apply_to_compound_groups() -> None:
    def dummy_fn(compound_group: List[str]) -> List[str]:
        # as a test: reverse the order of the compounds
        return list(reversed(compound_group))

    assert apply_to_compound_groups(
        ReactionEquation(["C", "N", "O"], ["P"], ["CNO", "C"]), dummy_fn
    ) == ReactionEquation(["O", "N", "C"], ["P"], ["C", "CNO"])


def test_rxn_standardization() -> None:
    initial_reaction_equation = ReactionEquation(
        reactants=["OCC", "O", "C"],
        agents=["CCO"],
        products=["C"],
    )

    # NB:
    # 1. CCO is merged into the reactants
    # 2. OCC is canonicalized to CCO
    # 3. CCO is now present two times and removed
    # 4. sorted alphabetically
    expected = ReactionEquation(
        reactants=["C", "CCO", "O"],
        agents=[],
        products=["C"],
    )

    assert rxn_standardization(initial_reaction_equation) == expected


def test_remove_precursors_from_products() -> None:
    reactants = ["COCO", "[Na+].[OH-]", "OCC"]
    agents = ["O", "C"]
    products = ["NCOC", "C", "COCO"]
    reaction = ReactionEquation(reactants, agents, products)

    reaction = remove_precursors_from_products(reaction)

    expected_products = ["NCOC"]
    expected_reaction = ReactionEquation(reactants, agents, expected_products)

    assert reaction == expected_reaction


def test_can_be_instantiated_from_any_iterator() -> None:
    a = {"C", "O"}
    b = (m for m in ["CO", "OC"])
    reaction_equation = ReactionEquation(a, [], b)

    assert reaction_equation.reactants == ["C", "O"] or reaction_equation.reactants == [
        "O",
        "C",
    ]
    assert reaction_equation.products == ["CO", "OC"]


def test_list_objects_are_not_shared() -> None:
    a = ["C", "O"]
    b = ["CO"]
    reaction_equation = ReactionEquation(a, [], b)

    # Modifying a list after using it to instantiate ReactionEquation should
    # have no impact on the ReactionEquation instance
    assert reaction_equation.reactants == a
    a.append("N")
    assert reaction_equation.reactants != a

    # Modifying the ReactionEquation compounds lists should have no impact on
    # the lists that were used to instantiate it.
    assert b == reaction_equation.products
    reaction_equation.products.append("N")
    assert b != reaction_equation.products
