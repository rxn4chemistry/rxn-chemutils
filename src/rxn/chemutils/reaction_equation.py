from functools import partial
from typing import (
    Callable,
    Generator,
    Iterable,
    Iterator,
    List,
    Optional,
    Type,
    TypeVar,
)

import attr
from rxn.utilities.containers import remove_duplicates

from .conversion import canonicalize_smiles, cleanup_smiles
from .exceptions import InvalidReactionSmiles
from .multicomponent_smiles import (
    list_to_multicomponent_smiles,
    multicomponent_smiles_to_list,
)

T = TypeVar("T", bound="ReactionEquation")


@attr.s(auto_attribs=True, init=False)
class ReactionEquation:
    """
    Defines a reaction equation, as given by the molecules involved in a reaction.

    Attributes:
        reactants: SMILES strings for compounds on the left of the reaction arrow.
        agents: SMILES strings for compounds above the reaction arrow. Are
            sometimes merged with the reactants.
        products: SMILES strings for compounds on the right of the reaction arrow.
    """

    reactants: List[str]
    agents: List[str]
    products: List[str]

    def __init__(
        self, reactants: Iterable[str], agents: Iterable[str], products: Iterable[str]
    ):
        """
        Overwrite init function in order to enable instantiation from any iterator and
        to force copying the lists.
        """
        self.__attrs_init__(  # type: ignore
            list(reactants), list(agents), list(products)
        )

    def __iter__(self) -> Iterator[List[str]]:
        """Helper function to simplify functionality acting on all three
        compound groups"""
        return (i for i in (self.reactants, self.agents, self.products))

    def iter_all_smiles(self) -> Generator[str, None, None]:
        """Helper function to iterate over all the SMILES in the reaction equation"""
        return (molecule for group in self for molecule in group)

    def to_string(self, fragment_bond: Optional[str] = None) -> str:
        """
        Convert a ReactionEquation to an "rxn" reaction SMILES.
        """

        smiles_groups = (
            list_to_multicomponent_smiles(group, fragment_bond) for group in self
        )
        return ">".join(smiles_groups)

    @classmethod
    def from_string(
        cls: Type[T], reaction_string: str, fragment_bond: Optional[str] = None
    ) -> T:
        """
        Convert a ReactionEquation from an "rxn" reaction SMILES.
        """

        groups = [
            multicomponent_smiles_to_list(smiles_group, fragment_bond=fragment_bond)
            for smiles_group in reaction_string.split(">")
        ]

        try:
            return cls(*groups)
        except TypeError as e:
            raise InvalidReactionSmiles(reaction_string) from e


def merge_reactants_and_agents(reaction: ReactionEquation) -> ReactionEquation:
    return ReactionEquation(
        reactants=reaction.reactants + reaction.agents,
        agents=[],
        products=reaction.products,
    )


def sort_compounds(reaction: ReactionEquation) -> ReactionEquation:
    """
    Reorder the compounds of each group in alphabetic order.
    """
    sorted_compound_groups = (sorted(group) for group in reaction)
    return ReactionEquation(*sorted_compound_groups)


def apply_to_compounds(
    reaction: ReactionEquation, fn: Callable[[str], str]
) -> ReactionEquation:
    """
    Apply a function to the individual compounds in a ReactionEquation.

    Args:
        reaction: reaction equation to apply the function to.
        fn: function to apply.

    Returns:
        New ReactionEquation instance after application of the function to the compounds.
    """
    updated_compound_groups = (
        [fn(compound) for compound in compound_group] for compound_group in reaction
    )
    return ReactionEquation(*updated_compound_groups)


def canonicalize_compounds(
    reaction: ReactionEquation, check_valence: bool = True
) -> ReactionEquation:
    """
    Canonicalize the molecules of a ReactionEquation.
    """
    canonicalize_fn = partial(canonicalize_smiles, check_valence=check_valence)
    return apply_to_compounds(reaction, canonicalize_fn)


def remove_duplicate_compounds(reaction: ReactionEquation) -> ReactionEquation:
    """
    Remove compounds that are duplicated in the same category
    """
    groups_without_duplicates = (remove_duplicates(group) for group in reaction)
    return ReactionEquation(*groups_without_duplicates)


def cleanup_compounds(reaction: ReactionEquation) -> ReactionEquation:
    """
    Basic cleanup of the compounds.
    """
    return apply_to_compounds(reaction, cleanup_smiles)


def rxn_standardization(reaction: ReactionEquation) -> ReactionEquation:
    """
    Apply the standard RXN postprocessing of reaction equations.

    Consists in the following
    1. merge reactants and agents
    2. canonicalize all the SMILES
    3. sort the compounds in each group
    4. remove the duplicates
    """
    return remove_duplicate_compounds(
        sort_compounds(canonicalize_compounds(merge_reactants_and_agents(reaction)))
    )


def remove_precursors_from_products(reaction: ReactionEquation) -> ReactionEquation:
    """
    Remove compounds in products that are also present in the reactants or reagents.
    """
    precursors = reaction.reactants + reaction.agents
    products_without_precursors = [
        product for product in reaction.products if product not in precursors
    ]
    return ReactionEquation(
        reactants=reaction.reactants,
        agents=reaction.agents,
        products=products_without_precursors,
    )


def has_repeated_molecules(reaction_equation: ReactionEquation) -> bool:
    all_molecules = list(reaction_equation.iter_all_smiles())
    return len(set(all_molecules)) < len(all_molecules)
