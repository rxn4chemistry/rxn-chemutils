import logging
import random
from typing import Callable, List

from .miscellaneous import apply_to_any_smiles, apply_to_smiles_groups

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class SmilesAugmenter:
    """
    Class to augment any kind of SMILES string with the help of randomization
    and shuffling.
    """

    def __init__(
        self,
        augmentation_fn: Callable[[str], str],
        augmentation_probability: float = 1.0,
        shuffle: bool = True,
        ignore_exceptions: bool = True,
    ):
        """
        Args:
            augmentation_fn: Function for augmenting the individual SMILES strings,
                such as the functions provided in smiles_randomization.py.
            augmentation_probability: Probability with which to augment individual
                SMILES strings.
            shuffle: Whether to shuffle the order of the compounds.
            ignore_exceptions: Whether to ignore the error (and return the
                original string) when an augmentation fails. If False, exceptions
                will be propagated.
        """
        self.augmentation_fn = augmentation_fn
        self.augmentation_probability = augmentation_probability
        self.shuffle = shuffle
        self.ignore_exceptions = ignore_exceptions

    def augment(self, smiles: str, number_augmentations: int) -> List[str]:
        """
        Augment one SMILES string (of any kind).

        Args:
            smiles: SMILES string to augment.
            number_augmentations: how many times to do the augmentation.
        """

        # augmentation of the individual compound SMILES
        augmented = [
            apply_to_any_smiles(
                smiles, self._augment_with_probability, force_multicomponent=True
            )
            for _ in range(number_augmentations)
        ]

        # shuffle the order of the compounds
        if self.shuffle:
            augmented = [
                apply_to_smiles_groups(s, SmilesAugmenter._shuffle) for s in augmented
            ]

        return augmented

    def _augment_with_probability(self, smiles: str) -> str:
        """Augmentat a SMILES, with the probability given by the member variable."""

        # Note: no need to call random.uniform if the augmentation probability is 1.0.
        if (
            self.augmentation_probability == 1.0
            or random.uniform(0, 1) <= self.augmentation_probability
        ):
            try:
                return self.augmentation_fn(smiles)
            except Exception as e:
                if self.ignore_exceptions:
                    logger.warning(f"Augmentation failed for {smiles}: {e}")
                    return smiles
                else:
                    raise

        # no augmentation
        return smiles

    @staticmethod
    def _shuffle(smiles_list: List[str]) -> List[str]:
        smiles_list = smiles_list.copy()
        random.shuffle(smiles_list)
        return smiles_list
