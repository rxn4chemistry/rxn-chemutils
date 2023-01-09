import random

import pytest
from rxn.utilities.basic import identity
from rxn.utilities.containers import all_identical

from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.smiles_augmenter import SmilesAugmenter
from rxn.chemutils.smiles_randomization import randomize_smiles_rotated

single_compound = "O=C(C)Oc1ccccc1C(=O)O"
salt_compound = "CC[NH+](CC)CC.CC(=O)[O-]"
multismiles = "O=C(C)Oc1ccccc1C(=O)O.Cl.[Na+]~[Cl-]"
rxn_smiles_1 = "CC(=O)Cl.NCC.[Na+]~[Cl-]>>CC(=O)NCC"
rxn_smiles_2 = "CC(=O)Cl.NCC.[Na+].[Cl-]>>CC(=O)NCC |f:2.3|"  # extended format
rxn_smiles_3 = "CC(=O)Cl.NCC.[Na+]~[Cl-]>>CC(=O)NCC.Cl"  # multiple products


def test_do_nothing() -> None:
    # If there is no SMILES randomization and no shuffling, the class should
    # just return the same strings back
    augmenter = SmilesAugmenter(augmentation_fn=identity, shuffle=False)

    test_strings = [
        single_compound,
        salt_compound,
        multismiles,
        rxn_smiles_1,
        rxn_smiles_2,
        rxn_smiles_3,
    ]
    for test_smiles in test_strings:
        assert augmenter.augment(test_smiles, 4) == 4 * [test_smiles]


def test_shuffling_only() -> None:
    augmenter = SmilesAugmenter(augmentation_fn=identity, shuffle=True)

    # On single compound: no change at all
    assert augmenter.augment(single_compound, 4) == 4 * [single_compound]

    # On salt compound: just reorders it -> 2 possible combinations
    assert len(set(augmenter.augment(salt_compound, 20))) == 2

    # On multismiles: fragments stay together; this leads to a total of 6 combinations
    assert len(set(augmenter.augment(multismiles, 60))) == 6

    # On reaction SMILES with two products: 6 combinations for precursors, 2 for products
    assert len(set(augmenter.augment(rxn_smiles_3, 150))) == 12


def test_smiles_augmentation_only() -> None:
    def dummy_augmentation(smiles: str) -> str:
        # Replace the SMILES strings by their lengths.
        return str(len(smiles))

    augmenter = SmilesAugmenter(augmentation_fn=dummy_augmentation, shuffle=False)

    query_and_expected = [
        (single_compound, "21"),
        (salt_compound, "13.10"),
        (multismiles, "21.2.11"),
        (rxn_smiles_1, "8.3.11>>9"),
        (rxn_smiles_2, "8.3.11>>9"),
        (rxn_smiles_3, "8.3.11>>9.2"),
    ]
    for query, expected in query_and_expected:
        # Note: no randomness, as the dummy augmentation is not random, and
        # there is no shuffling
        assert augmenter.augment(query, 4) == 4 * [expected]


def test_smiles_augmentation_only_with_probability() -> None:
    def dummy_augmentation(smiles: str) -> str:
        # Replace the SMILES strings by their lengths.
        return str(len(smiles))

    augmenter = SmilesAugmenter(
        augmentation_fn=dummy_augmentation, augmentation_probability=0.5, shuffle=False
    )

    # On single compound: either augmented or not
    assert set(augmenter.augment(single_compound, 10)) == {"21", single_compound}

    # To illustrate the behavior on strings with multiple SMILES, we just test it on
    # the multismiles input. There will be a mix of replaced and not replaced
    assert set(augmenter.augment(multismiles, 50)) == {
        "O=C(C)Oc1ccccc1C(=O)O.Cl.[Na+]~[Cl-]",
        "21.Cl.[Na+]~[Cl-]",
        "O=C(C)Oc1ccccc1C(=O)O.2.[Na+]~[Cl-]",
        "O=C(C)Oc1ccccc1C(=O)O.Cl.11",
        "21.2.[Na+]~[Cl-]",
        "21.Cl.11",
        "O=C(C)Oc1ccccc1C(=O)O.2.11",
        "21.2.11",
    }

    # No augmentation at all if the probability is zero
    augmenter = SmilesAugmenter(
        augmentation_fn=dummy_augmentation, augmentation_probability=0.0, shuffle=False
    )
    assert set(augmenter.augment(multismiles, 10)) == {multismiles}


def test_mix_augmentation_and_shuffling() -> None:
    augmenter = SmilesAugmenter(
        augmentation_fn=randomize_smiles_rotated,
        augmentation_probability=0.5,
        shuffle=True,
    )

    random.seed(42)
    assert augmenter.augment(rxn_smiles_3, 6) == [
        "NCC.CC(=O)Cl.[Cl-]~[Na+]>>Cl.CC(=O)NCC",
        "NCC.ClC(=O)C.[Na+]~[Cl-]>>Cl.N(CC)C(C)=O",
        "C(=O)(Cl)C.[Na+]~[Cl-].CCN>>CC(=O)NCC.Cl",
        "NCC.CC(=O)Cl.[Na+]~[Cl-]>>CC(=O)NCC.Cl",
        "NCC.[Na+]~[Cl-].CC(=O)Cl>>Cl.CC(=O)NCC",
        "CC(=O)Cl.[Na+]~[Cl-].C(C)N>>CC(=O)NCC.Cl",
    ]


def test_reproducibility() -> None:
    # illustrate the reproducibility at the example of a reaction SMILES
    augmenter = SmilesAugmenter(
        augmentation_fn=randomize_smiles_rotated,
        augmentation_probability=0.5,
        shuffle=True,
    )

    # When resetting the seed, we should always get exactly the same results
    results = []
    for _ in range(10):
        random.seed(42)
        results.append(augmenter.augment(rxn_smiles_3, 5))
    assert all_identical(results)

    # sampling one more time without resetting the seed -> results change
    results.append(augmenter.augment(rxn_smiles_3, 5))
    assert not all_identical(results)


def test_augmentation_errors() -> None:
    augmenter = SmilesAugmenter(
        augmentation_fn=randomize_smiles_rotated, ignore_exceptions=True
    )

    invalid_smiles = "thisisinvalid"

    # When errors are ignored: returns the original input
    assert augmenter.augment(invalid_smiles, 1) == [invalid_smiles]

    # When errors are not ignored: raises an exception
    augmenter.ignore_exceptions = False
    with pytest.raises(InvalidSmiles):
        _ = augmenter.augment(invalid_smiles, 1)
