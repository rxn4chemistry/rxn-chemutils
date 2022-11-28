from rxn.utilities.basic import identity

from rxn.chemutils.extended_reaction_smiles import to_extended_reaction_smiles
from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles
from rxn.chemutils.smiles_augmenter import SmilesAugmenter

single_compound = "O=C(C)Oc1ccccc1C(=O)O"
salt_compound = "CC[NH+](CC)CC.CC(=O)[O-]"
multismiles = "O=C(C)Oc1ccccc1C(=O)O.Cl.[Na+]~[Cl-]"
rxn_smiles_1 = "CC(=O)Cl.NCC.[Na+]~[Cl-]>>CC(=O)NCC"
rxn_smiles_2 = "CC(=O)Cl.NCC.[Na+].[Cl-]>>CC(=O)NCC |f:2.3|"  # extended format
rxn_smiles_3 = "CC(=O)Cl.NCC.[Na+]~[Cl-]>>CC(=O)NCC.Cl"  # multiple products


def test_do_nothing() -> None:
    # If there is no SMILES randomization and no shuffling, the class should
    # just return the same strings back
    augmenter = SmilesAugmenter(augmentation_fn=identity, shuffle_order=False)

    test_strings = [
        single_compound,
        salt_compound,
        multismiles,
        rxn_smiles_1,
        rxn_smiles_2,
        rxn_smiles_3,
    ]
    for test_smiles in test_strings:
        assert augmenter.augment_one(test_smiles, 4) == 4 * [test_smiles]


def test_shuffling_only() -> None:
    augmenter = SmilesAugmenter(augmentation_fn=identity, shuffle_order=True)


def test_exceptions() -> None:
    pass
