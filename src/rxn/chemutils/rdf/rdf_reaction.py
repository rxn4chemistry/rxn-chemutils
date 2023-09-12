from typing import Dict, List

from attr import define


@define
class RdfReaction:
    """
    Custom class wrapping an RDF reaction.

    Attributes:
        reactants: List of reactants, given in MDL format.
        reagents: List of reagents, given in MDL format.
        products: List of products, given in MDL format.
        meta: meta information.
        reaction_index: index of the reaction, coming from the "$RFMT $RIREG" line.
    """

    reactants: List[str]
    reagents: List[str]
    products: List[str]
    meta: Dict[str, str]
    reaction_index: int

    def get_mol_structure(self, category: str, index: int) -> str:
        """
        Get a structure given in the reaction metadata, such as in
        RXN:CATALYST(2):MOL(1):MOLSTRUCTURE

        Args:
            category: 'catalyst' for the example above.
            index: 2 for the example above.

        Returns:
            The MolBlock for the desired structure.
        """
        key = f"RXN:{category.upper()}({index}):MOL(1):MOLSTRUCTURE"
        return self.meta[key]

    def additional_reagents_from_meta(self) -> List[str]:
        """
        Get the MDL strings for the additional reagents from the meta information.
        """
        return [
            value for key, value in self.meta.items() if key.endswith("MOLSTRUCTURE")
        ]

    def all_reagents(self) -> List[str]:
        """
        Get the MDL strings for all the reagents.
        """

        reagent_mdls = self.reagents.copy()
        reagent_mdls.extend(self.additional_reagents_from_meta())
        return reagent_mdls
