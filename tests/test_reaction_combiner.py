import pytest

from rxn.chemutils.reaction_combiner import ReactionCombiner
from rxn.chemutils.reaction_smiles import ReactionFormat


def test_reaction_combiner_on_precursors_and_products() -> None:
    combiner = ReactionCombiner()

    precursors = ["CC.O", "CCC.O"]
    products = ["CCO", "CCCO"]
    expected = ["CC.O>>CCO", "CCC.O>>CCCO"]

    assert list(combiner.combine(precursors, products)) == expected


def test_reaction_combiner_on_fragment_reactions() -> None:
    combiner = ReactionCombiner()

    fragment_1 = ["CC.O>>", "CCC>>CCCO"]
    fragment_2 = [">>CCO", "O.N>>"]
    expected = ["CC.O>>CCO", "CCC.O.N>>CCCO"]

    assert list(combiner.combine(fragment_1, fragment_2)) == expected


def test_reaction_combiner_with_tokenized_input() -> None:
    combiner = ReactionCombiner()

    fragment_1 = ["C C . O", "C C C >> C C C O"]
    fragment_2 = ["C C O", "O . N >>"]
    expected = ["CC.O>>CCO", "CCC.O.N>>CCCO"]

    assert list(combiner.combine(fragment_1, fragment_2)) == expected


def test_multiple_precursors_per_product() -> None:
    combiner = ReactionCombiner()

    # Three sets of precursors for one product
    # Note: the same would work for two sets of fragment reactions.
    precursors = ["CC.O", "CC.O.N", "CC.O.P", "CCC.O", "CCC.O.N", "CCC.O.P"]
    products = ["CCO", "CCCO"]
    expected = [
        "CC.O>>CCO",
        "CC.O.N>>CCO",
        "CC.O.P>>CCO",
        "CCC.O>>CCCO",
        "CCC.O.N>>CCCO",
        "CCC.O.P>>CCCO",
    ]

    assert list(combiner.combine(precursors, products)) == expected


def test_multiple_products_per_precursors() -> None:
    combiner = ReactionCombiner()

    # Two products for each set of precursors
    # Note: the same would work for two sets of fragment reactions.
    precursors = ["CC.O", "CCC.O"]
    products = ["CCO", "OCCO", "CCCO", "OCCCO"]
    expected = [
        "CC.O>>CCO",
        "CC.O>>OCCO",
        "CCC.O>>CCCO",
        "CCC.O>>OCCCO",
    ]

    assert list(combiner.combine(precursors, products)) == expected


def test_incompatible_number_of_precursors_and_products() -> None:
    combiner = ReactionCombiner()

    # Three precursors and two products -> unclear how to combine
    precursors = ["CC.O", "CCC.O", "CCCC.O"]
    products = ["CCO", "CCCO"]

    with pytest.raises(ValueError):
        _ = list(combiner.combine(precursors, products))


def test_different_output_formats() -> None:
    fragment_1 = ["CC~O", "CCC>>CCCO"]
    fragment_2 = ["CCO", "N~O>>"]

    assert list(
        ReactionCombiner(reaction_format=ReactionFormat.STANDARD).combine(
            fragment_1, fragment_2
        )
    ) == ["CC.O>>CCO", "CCC.N.O>>CCCO"]
    assert list(
        ReactionCombiner(reaction_format=ReactionFormat.STANDARD_WITH_TILDE).combine(
            fragment_1, fragment_2
        )
    ) == ["CC~O>>CCO", "CCC.N~O>>CCCO"]
    assert list(
        ReactionCombiner(reaction_format=ReactionFormat.EXTENDED).combine(
            fragment_1, fragment_2
        )
    ) == ["CC.O>>CCO |f:0.1|", "CCC.N.O>>CCCO |f:1.2|"]


def test_standardization() -> None:
    precursors = ["C(C).O", "O.C(CC)"]
    products = ["CCO", "CCCO"]

    assert list(ReactionCombiner(standardize=True).combine(precursors, products)) == [
        "CC.O>>CCO",
        "CCC.O>>CCCO",
    ]
    assert list(ReactionCombiner(standardize=False).combine(precursors, products)) == [
        "C(C).O>>CCO",
        "O.C(CC)>>CCCO",
    ]


def test_invalid_reactions() -> None:
    fallback_string = "C>>C"
    combiner = ReactionCombiner(standardize=True, fallback_reaction=fallback_string)

    precursors = ["CC.O", "CCC.O", "CC>CC>O"]
    products = ["CCo", "CCCO", "CCCCO"]
    expected = [fallback_string, "CCC.O>>CCCO", fallback_string]

    assert list(combiner.combine(precursors, products)) == expected


def test_combine_iterators() -> None:
    combiner = ReactionCombiner()

    precursors = iter(["CC.O", "CCC.O"])
    products = iter(["CCO", "CCCO"])
    expected = ["CC.O>>CCO", "CCC.O>>CCCO"]

    assert list(combiner.combine_iterators(precursors, products)) == expected


def test_combine_iterators_with_multiplier() -> None:
    combiner = ReactionCombiner()

    precursors = iter(["CC.O", "CC.O.N", "CCC.O", "CCC.O.N"])
    products = iter(["CCO", "CCCO"])
    expected = ["CC.O>>CCO", "CC.O.N>>CCO", "CCC.O>>CCCO", "CCC.O.N>>CCCO"]

    results = combiner.combine_iterators(precursors, products, 1, 2)
    assert list(results) == expected
