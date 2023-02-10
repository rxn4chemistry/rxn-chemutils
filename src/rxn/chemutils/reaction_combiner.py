from itertools import chain, repeat
from typing import Iterator, List, Tuple

from rxn.utilities.misc import get_multipliers

from .reaction_equation import canonicalize_compounds, sort_compounds
from .reaction_smiles import (
    ReactionFormat,
    parse_any_reaction_smiles,
    to_reaction_smiles,
)
from .tokenization import detokenize_smiles


class ReactionCombiner:
    """
    Class to combine sets of precursors with sets of products, for instance to
    merge src and prediction files.
    """

    def __init__(
        self,
        standardize: bool = False,
        reaction_format: ReactionFormat = ReactionFormat.STANDARD_WITH_TILDE,
        invalid_reaction: str = ">>",
    ):
        """
        Args:
            standardize: whether to standardize (i.e. canonicalize and reorder) the reaction SMILES.
            reaction_format: which format should be used for the reaction SMILES.
            invalid_reaction: text / reaction to produce when a reaction is invalid.
        """
        self.standardize = standardize
        self.reaction_format = reaction_format
        self.invalid_reaction = invalid_reaction

    def combine(self, precursors: List[str], products: List[str]) -> Iterator[str]:
        """
        Args:
            precursors: List of sets of precursors strings (such as "CC.O.[Na+]~[Cl-]").
            products: List of sets of product(s) strings.

        Returns:
            Iterator over the resulting reaction SMILES.
        """
        precursor_multiplier, product_multiplier = self._get_multipliers(
            precursors, products
        )

        # repeat itemwise the elements: https://stackoverflow.com/a/45799308
        precursor_iterator = chain.from_iterable(
            zip(*repeat(precursors, precursor_multiplier))
        )
        product_iterator = chain.from_iterable(
            zip(*repeat(products, product_multiplier))
        )

        for precursor_string, product_string in zip(
            precursor_iterator, product_iterator
        ):
            yield self._to_reaction_smiles(precursor_string, product_string)

    def _to_reaction_smiles(self, precursor_string: str, product_string: str) -> str:
        try:
            return self._try_to_reaction_smiles(precursor_string, product_string)
        except Exception:
            return self.invalid_reaction

    def _try_to_reaction_smiles(
        self, precursor_string: str, product_string: str
    ) -> str:
        # 1) get the initial reaction SMILES
        reaction_smiles = precursor_string + ">>" + product_string
        reaction_smiles = detokenize_smiles(reaction_smiles)

        # 2) standardize if necessary
        reaction_equation = parse_any_reaction_smiles(reaction_smiles)
        if self.standardize:
            reaction_equation = sort_compounds(
                canonicalize_compounds(reaction_equation)
            )

        return to_reaction_smiles(
            reaction_equation, reaction_format=self.reaction_format
        )

    def _get_multipliers(
        self, precursors: List[str], products: List[str]
    ) -> Tuple[int, int]:
        """Get the multipliers to use when iterating through the precursors and products.

        Raises:
            ValueError: when one is not exactly a multiple of the other.

        Returns:
            Tuple: precursor multiplier, product multiplier
        """
        a = len(precursors)
        b = len(products)

        m_a, m_b = get_multipliers(a, b)

        # Fail if one is not exactly a multiple of the other
        if 1 not in {m_a, m_b}:
            raise ValueError(
                "The number of precursor sets and the number of products "
                f"are not an exact multiple of each other ({a} and {b})"
            )

        return m_a, m_b
