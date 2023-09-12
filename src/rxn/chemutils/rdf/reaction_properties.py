from typing import Any, Dict, Iterator, Optional, Tuple

import attr

from .property_parser import parse_properties


@attr.s(auto_attribs=True)
class Compound:
    """Wrapper for the compound dicts, which always contain a MOLSTRUCTURE,
    usually a SYMBOL, and sometimes a MOL_CAPTION."""

    mol_structure: str
    symbol: Optional[str]
    mol_caption: Optional[str]
    category: Optional[str]

    def get_name(self) -> Optional[str]:
        """Try to get a name by first looking at the SYMBOL and then at the MOLCAPTION."""
        if self.symbol is not None:
            return self.symbol
        return self.mol_caption

    def to_dict(self) -> Dict[str, str]:
        """Convert back to a dictionary, as in the RdfReaction properties."""
        d = {
            "MOLSTRUCTURE": self.mol_structure,
        }
        if self.symbol is not None:
            d["SYMBOL"] = self.symbol
        if self.mol_caption is not None:
            d["MOL_CAPTION"] = self.mol_caption
        return d

    @classmethod
    def from_dict(
        cls, compound_dict: Dict[str, str], category: Optional[str] = None
    ) -> "Compound":
        if category == "":
            category = None

        return cls(
            mol_structure=compound_dict["MOLSTRUCTURE"],
            symbol=compound_dict.get("SYMBOL", None),
            mol_caption=compound_dict.get("MOL_CAPTION", None),
            category=category,
        )


class ReactionProperties:
    """
    Class compiling the reaction properties given for an RDF reaction.
    """

    def __init__(self, meta: Dict[str, str]):
        pruned_meta = self._prune_unnecessary_lists(meta)
        parsed_properties = parse_properties(pruned_meta)
        try:
            self.properties: Dict[str, Any] = parsed_properties["RXN"]
        except KeyError:
            self.properties = parsed_properties

    def get_compound_dicts(self) -> Iterator[Dict[str, str]]:
        """NB: this function allows to iterate through the compounds in the
        property dict and make edits there (see ReactionUpdater)."""
        yield from find_compounds(self.properties)

    def get_compounds(self) -> Iterator[Compound]:
        """Get the compounds involved in a reaction, without the reactants and products."""
        yield from (
            Compound.from_dict(cpd_dict, category)
            for cpd_dict, category in find_compounds_with_category(self.properties)
        )

    def _prune_unnecessary_lists(self, meta: Dict[str, str]) -> Dict[str, str]:
        """
        In the meta dict, often there will be unnecessary list-like keys.
        In particular, "MOL(1):" and "SYMBOL(1):" can be removed as they are
        never necessary.
        """

        def simplify_key(key: str) -> str:
            updated_key = key.replace("MOL(1):", "").replace("SYMBOL(1):", "")
            if "MOL(" in updated_key or "SYMBOL(" in updated_key:
                raise ValueError(f"MOL or SYMBOL list with index other than 1: {key}")
            return updated_key

        return {simplify_key(key): value for key, value in meta.items()}


def find_compounds(property_dict: Dict[str, Any]) -> Iterator[Dict[str, str]]:
    """Iterate through a property dictionary to find compound dicts."""
    yield from (
        compound_dict
        for compound_dict, _ in find_compounds_with_category(property_dict)
    )


def find_compounds_with_category(
    property_dict: Dict[str, Any]
) -> Iterator[Tuple[Dict[str, str], str]]:
    """Iterate through a property dictionary to find compound dicts, with the
    associated category ("SOLVENT", "CATALYST", etc.)."""
    yield from _find_compounds_with_prefix("", property_dict)


def _find_compounds_with_prefix(
    prefix: str, property_dict: Dict[str, Any]
) -> Iterator[Tuple[Dict[str, str], str]]:
    if "MOLSTRUCTURE" in property_dict:
        yield property_dict, prefix

    for key, value in property_dict.items():
        if isinstance(value, dict):
            yield from _find_compounds_with_prefix("", value)

        # If the item is a list, we take the key as a prefix - this may be
        # something such as "SOLVENT", "CATALYST", etc.
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    yield from _find_compounds_with_prefix(key, item)
