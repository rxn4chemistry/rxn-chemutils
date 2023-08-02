from itertools import chain, repeat, zip_longest
from typing import Iterable, Iterator, Sequence, Tuple

from rxn.utilities.misc import get_multipliers

from .miscellaneous import merge_reactions
from .reaction_equation import ReactionEquation, canonicalize_compounds, sort_compounds
from .reaction_smiles import (
    ReactionFormat,
    parse_any_reaction_smiles,
    to_reaction_smiles,
)
from .tokenization import detokenize_smiles


class ReactionCombiner:
    """
    Class to combine sets of precursors with sets of products, or sets of partial
    reactions with other sets of partial reactions.

    This class is typically useful when one needs to produce the full reaction
    SMILES starting from multiple files, such as A) one file for the precursors
    and one for the products, or B) two files containing each one part of a
    chemical equation.

    This class is particularly useful when the said files have different sizes,
    which can be the case when multiple predictions are made for each line of one
    of these files.
    """

    def __init__(
        self,
        standardize: bool = False,
        reaction_format: ReactionFormat = ReactionFormat.STANDARD_WITH_TILDE,
        fallback_reaction: str = ">>",
    ):
        """
        Args:
            standardize: whether to standardize (i.e. canonicalize and reorder) the reaction SMILES.
            reaction_format: which format should be used for the reaction SMILES.
            fallback_reaction: text / reaction to produce when a reaction is invalid.
        """
        self.standardize = standardize
        self.reaction_format = reaction_format
        self.fallback_reaction = fallback_reaction

    def combine(
        self, fragments_1: Sequence[str], fragments_2: Sequence[str]
    ) -> Iterator[str]:
        """
        See docstring of function ``combine_sequences``.
        """
        yield from self.combine_sequences(fragments_1, fragments_2)

    def combine_sequences(
        self, fragments_1: Sequence[str], fragments_2: Sequence[str]
    ) -> Iterator[str]:
        """
        Combine the two sequences of fragments into an iterator of reactions.

        Args:
            fragments_1: Sequence of sets of precursors strings (such as "CC.O.[Na+]~[Cl-]"),
                or list of partial reactions.
            fragments_2: Sequence of sets of product(s) strings, or list of partial
                reactions.

        Returns:
            Iterator over the resulting reaction SMILES.
        """
        fragments_1_multiplier, fragments_2_multiplier = self._get_multipliers(
            fragments_1, fragments_2
        )
        yield from self.combine_iterators(
            fragments_1=fragments_1,
            fragments_2=fragments_2,
            fragments_1_multiplier=fragments_1_multiplier,
            fragments_2_multiplier=fragments_2_multiplier,
        )

    def combine_iterators(
        self,
        fragments_1: Iterable[str],
        fragments_2: Iterable[str],
        fragments_1_multiplier: int = 1,
        fragments_2_multiplier: int = 1,
    ) -> Iterator[str]:
        """
        Combine the two iterators of fragments into an iterator of reactions.

        Args:
            fragments_1: Sequence of sets of precursors strings (such as "CC.O.[Na+]~[Cl-]"),
                or list of partial reactions.
            fragments_2: Sequence of sets of product(s) strings, or list of partial
                reactions.
            fragments_1_multiplier: how many times to duplicate the fragments_1.
            fragments_2_multiplier: how many times to duplicate the fragments_2.

        Raises:
            RuntimeError: if one of the iterators isn't fully consumed.
            ValueError: when one is not exactly a multiple of the other.

        Returns:
            Iterator over the resulting reaction SMILES.
        """

        self._validate_multipliers(fragments_1_multiplier, fragments_2_multiplier)

        # repeat itemwise the elements: https://stackoverflow.com/a/45799320
        fragment_1_iterator = chain.from_iterable(
            (repeat(e, fragments_1_multiplier) for e in fragments_1)
        )
        fragment_2_iterator = chain.from_iterable(
            (repeat(e, fragments_2_multiplier) for e in fragments_2)
        )

        for fragment_1, fragment_2 in zip_longest(
            fragment_1_iterator, fragment_2_iterator
        ):
            if fragment_1 is None or fragment_2 is None:
                raise RuntimeError("Mismatch in expected iterator length")
            yield self._to_reaction_smiles(fragment_1, fragment_2)

    def _to_reaction_smiles(self, fragment_1: str, fragment_2: str) -> str:
        try:
            return self._try_to_reaction_smiles(fragment_1, fragment_2)
        except Exception:
            return self.fallback_reaction

    def _try_to_reaction_smiles(self, fragment_1: str, fragment_2: str) -> str:
        # 1) get the initial reaction SMILES
        reaction_equation = self._to_raw_reaction(fragment_1, fragment_2)

        # 2) standardize if necessary
        if self.standardize:
            reaction_equation = sort_compounds(
                canonicalize_compounds(reaction_equation)
            )

        return to_reaction_smiles(
            reaction_equation, reaction_format=self.reaction_format
        )

    def _to_raw_reaction(self, fragment_1: str, fragment_2: str) -> ReactionEquation:
        """Get a ReactionEquation from the two strings."""
        fragment_1 = detokenize_smiles(fragment_1)
        fragment_2 = detokenize_smiles(fragment_2)

        fragment_1_is_reaction = ">" in fragment_1
        fragment_2_is_reaction = ">" in fragment_2

        # Case A: both are given in the reaction format
        if fragment_1_is_reaction and fragment_2_is_reaction:
            reaction_1 = parse_any_reaction_smiles(fragment_1)
            reaction_2 = parse_any_reaction_smiles(fragment_2)
            return merge_reactions(reaction_1, reaction_2)

        # Case A: fragment_1 represents the precursor(s), fragment_2 the product(s)
        if not fragment_1_is_reaction and not fragment_2_is_reaction:
            reaction_smiles = fragment_1 + ">>" + fragment_2
            return parse_any_reaction_smiles(reaction_smiles)

        raise ValueError(
            f'Cannot determine how to combine "{fragment_1}" and "{fragment_2}"'
        )

    def _get_multipliers(
        self, fragments_1: Sequence[str], fragments_2: Sequence[str]
    ) -> Tuple[int, int]:
        """Get the multipliers to use when iterating through the respective fragments.

        Raises:
            ValueError: when one is not exactly a multiple of the other.

        Returns:
            Tuple: fragments_1 multiplier, fragments_2 multiplier
        """
        a = len(fragments_1)
        b = len(fragments_2)

        m_a, m_b = get_multipliers(a, b)

        return m_a, m_b

    def _validate_multipliers(self, multiplier_1: int, multiplier_2: int) -> None:
        """
        Make sure that the given multipliers can be used with the reaction combiner.

        Raises:
            ValueError: when one is not exactly a multiple of the other.
        """
        # Fail if one is not exactly a multiple of the other
        if 1 not in {multiplier_1, multiplier_2}:
            raise ValueError(
                "The number of fragments of reactions are not an exact multiple of "
                f"each other: the multipliers are {multiplier_1} and {multiplier_2}."
            )
