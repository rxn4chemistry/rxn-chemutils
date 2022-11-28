import random
from typing import Callable, Iterable, Iterator, List

from .miscellaneous import apply_to_any_smiles
from .reaction_equation import ReactionEquation
from .reaction_smiles import parse_any_reaction_smiles
from .smiles_randomization import (
    randomize_smiles_restricted,
    randomize_smiles_rotated,
    randomize_smiles_unrestricted,
)


class SmilesAugmenter:
    """
    Class to augment any kind of SMILES string with the help of randomization
    and shuffling.
    """

    def __init__(
        self,
        augmentation_fn: Callable[[str], str],
        augmentation_probability: float = 1.0,
        shuffle_order: bool = True,
    ):
        """
        Args:
            augmentation_fn: Function for augmenting the individual SMILES strings,
                such as the functions provided in smiles_randomization.py.
            augmentation_probability: Probability with which to augment individual
                SMILES strings.
            shuffle_order: Whether to shuffle the order of the compounds.
        """
        self.augmentation_fn = augmentation_fn
        self.augmentation_probability = augmentation_probability
        self.shuffle_order = shuffle_order

    def augment_one(self, smiles: str, number_augmentations: int) -> List[str]:
        """
        Augment one SMILES string (of any kind).

        Args:
            smiles: SMILES string to augment.
            number_augmentations: how many times to do the augmentation.
        """
        augmented = [
            apply_to_any_smiles(smiles, self.augmentation_fn)
            for _ in range(number_augmentations)
        ]
        return augmented

    # def augment_many(
    #     self, smiles_strings: Iterable[str], number_augmentations: int
    # ) -> Iterator[List[str]]:
    #     """
    #     Augment many SMILES strings (of any kind).
    #
    #     Args:
    #         smiles_strings: iterable over SMILES strings.
    #         number_augmentations: how many times to augment each SMILES string.
    #
    #     Returns:
    #         Iterator over the augmentations
    #     """
    #     for smiles in smiles_strings:
    #         yield self.augment_one(smiles, number_augmentations)


randomization_functions = [
    randomize_smiles_rotated,
    randomize_smiles_unrestricted,
    randomize_smiles_restricted,
]

fn = randomization_functions[0]

random.seed(42)

rxn_smiles = "CC(C)c1ccc(C(=O)CCCCl)cc1.ClC(Cl)(Cl)Cl.O=C(OOC(=O)c1ccccc1)c1ccccc1.O=C1CCC(=O)N1Br>>CC(C)(Br)c1ccc(C(=O)CCCCl)cc1"
print(rxn_smiles)

reaction: ReactionEquation = parse_any_reaction_smiles(rxn_smiles)


for _ in range(10):
    reactants = [fn(smiles) for smiles in reaction.reactants]
    products = [fn(smiles) for smiles in reaction.products]

    random.shuffle(reactants)
    random.shuffle(products)

    augmented_reaction = ReactionEquation(reactants, [], products)
    print(augmented_reaction.to_string(fragment_bond="~"))
