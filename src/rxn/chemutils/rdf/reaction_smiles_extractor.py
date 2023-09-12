from typing import Iterable, List, Optional

from rxn.chemutils.reaction_equation import ReactionEquation

from ..conversion import mdl_to_smiles
from .rdf_reaction import RdfReaction
from .reaction_properties import ReactionProperties


class ReactionSmilesExtractor:
    """Extract reaction SMILES from RdfReaction instances."""

    def __init__(self, fragment_bond: Optional[str] = None, sanitize: bool = True):
        self.fragment_bond = fragment_bond
        self.sanitize = sanitize

    def to_reaction_smiles(self, reaction: RdfReaction) -> str:
        """Extract the reaction SMILES.

        Args:
            reaction: RdfReaction to extract the reaction equation from.
        """
        return self.to_string(self.to_reaction_equation(reaction))

    def to_reaction_equation(self, reaction: RdfReaction) -> ReactionEquation:
        """Extract the reaction equation.

        Args:
            reaction: RdfReaction to extract the reaction equation from.
        """

        # NB: reaction.reagents should actually always be an empty list
        reagents = reaction.reagents.copy()
        reagents.extend(
            c.mol_structure for c in ReactionProperties(reaction.meta).get_compounds()
        )

        return ReactionEquation(
            reactants=self._to_smiles_group(reaction.reactants),
            agents=self._to_smiles_group(reagents),
            products=self._to_smiles_group(reaction.products),
        )

    def to_string(self, reaction_equation: ReactionEquation) -> str:
        """Convert a ReactionEquation to a reaction SMILES.

        Useful function as it automatically uses the correct fragment bond."""
        return reaction_equation.to_string(fragment_bond=self.fragment_bond)

    def _to_smiles(self, mdl: str) -> str:
        return mdl_to_smiles(mdl, sanitize=self.sanitize, canonicalize=True)

    def _to_smiles_group(self, mdl_iterable: Iterable[str]) -> List[str]:
        return [self._to_smiles(m) for m in mdl_iterable]
