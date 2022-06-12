"""
Functionality to handle different kinds of reaction SMILES with the same functions.

In a separate file than miscellaneous.py or conversion.py in order to avoid
cyclic dependencies.
"""
from enum import auto

from rxn.utilities.types import RxnEnum

from .extended_reaction_smiles import (
    parse_extended_reaction_smiles,
    to_extended_reaction_smiles,
)
from .reaction_equation import ReactionEquation


class ReactionFormat(RxnEnum):
    """
    Existing reaction SMILES formats.

    Attributes:
        EXTENDED: extended reaction SMILES with fragment info,
            f.i. "|f:0.2,3.4.5|".
        STANDARD: standard reaction SMILES.
        STANDARD_WITH_TILDE: standard reaction SMILES, where fragments are
            indicated with tilde symbols, "~".
    """

    EXTENDED = auto()
    STANDARD = auto()
    STANDARD_WITH_TILDE = auto()


def determine_format(reaction_smiles: str) -> ReactionFormat:
    """
    Determine the format of a reaction SMILES.
    """
    if " |f:" in reaction_smiles:
        return ReactionFormat.EXTENDED

    if "~" in reaction_smiles:
        return ReactionFormat.STANDARD_WITH_TILDE

    return ReactionFormat.STANDARD


def parse_any_reaction_smiles(smiles: str) -> ReactionEquation:
    """
    Parse a reaction SMILES in any format (will be determined automatically).
    """
    return parse_reaction_smiles(smiles, reaction_format=determine_format(smiles))


def parse_reaction_smiles(
    smiles: str, reaction_format: ReactionFormat
) -> ReactionEquation:
    """
    Parse the reaction SMILES in a given format.
    """
    if reaction_format is ReactionFormat.EXTENDED:
        return parse_extended_reaction_smiles(smiles, remove_atom_maps=False)

    if reaction_format is ReactionFormat.STANDARD:
        return ReactionEquation.from_string(smiles)

    if reaction_format is ReactionFormat.STANDARD_WITH_TILDE:
        return ReactionEquation.from_string(smiles, fragment_bond="~")

    raise ValueError(f"Unsupported reaction format: {reaction_format}")


def to_reaction_smiles(
    reaction_equation: ReactionEquation, reaction_format: ReactionFormat
) -> str:
    """
    Convert a reaction equation into a reaction SMILES of the specified format.
    """
    if reaction_format is ReactionFormat.EXTENDED:
        return to_extended_reaction_smiles(reaction_equation)

    if reaction_format is ReactionFormat.STANDARD:
        return reaction_equation.to_string()

    if reaction_format is ReactionFormat.STANDARD_WITH_TILDE:
        return reaction_equation.to_string(fragment_bond="~")

    raise ValueError(f"Unsupported reaction format: {reaction_format}")
