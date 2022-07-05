import pytest

from rxn.chemutils.reaction_equation import ReactionEquation
from rxn.chemutils.reaction_smiles import (
    ReactionFormat,
    determine_format,
    parse_any_reaction_smiles,
    parse_reaction_smiles,
    to_reaction_smiles,
)


def test_determine_format():
    reaction_smiles_with_other_extended_information = (
        "[CH3:1][C:2](=[O:3])[NH:4][C@H:5]([C:6](=[O:7])[OH:8])[C:9]([CH3:10])"
        "([CH3:11])[SH:12].[O:13]=[N:14]O[Na]>O.Cl.CO>[CH3:1][C:2](=[O:3])[NH:4]"
        "[C@H:5]([C:6](=[O:7])[OH:8])[C:9]([CH3:10])([CH3:11])[S:12][N:14]=[O:13]"
        " |c:5,7,H:3.2,&1:4,24,f:3.4|"
    )

    assert determine_format("CC.O>>CCO") is ReactionFormat.STANDARD
    assert (
        determine_format("CC.O.[Na+]~[Cl-]>>CCO") is ReactionFormat.STANDARD_WITH_TILDE
    )
    assert determine_format("CC.O.[Na+].[Cl-]>>CCO |f:2.3|") is ReactionFormat.EXTENDED
    assert (
        determine_format(reaction_smiles_with_other_extended_information)
        is ReactionFormat.EXTENDED
    )


def test_parse_any_reaction_smiles():
    assert parse_any_reaction_smiles("CC.O>>CCO") == ReactionEquation(
        ["CC", "O"], [], ["CCO"]
    )
    assert parse_any_reaction_smiles("CC.O.[Na+]~[Cl-]>>CCO") == ReactionEquation(
        ["CC", "O", "[Na+].[Cl-]"], [], ["CCO"]
    )
    assert parse_any_reaction_smiles(
        "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"
    ) == ReactionEquation(["CC", "O", "[Na+].[Cl-]"], [], ["CCO"])


def test_parse_any_reaction_smiles_with_mapping():
    # Molecules with mapping
    reactant_1 = "[CH3:1][CH3:2]"
    reactant_2 = "[OH2:3]"
    product = "[CH3:1][CH2:2][OH:3]"

    # In the extracted equation, mapping is still present.
    expected = ReactionEquation([reactant_1, reactant_2, "[Na+].[Cl-]"], [], [product])

    assert (
        parse_any_reaction_smiles(f"{reactant_1}.{reactant_2}.[Na+]~[Cl-]>>{product}")
        == expected
    )
    assert (
        parse_any_reaction_smiles(
            f"{reactant_1}.{reactant_2}.[Na+].[Cl-]>>{product} |f:2.3|"
        )
        == expected
    )


def test_parse_reaction_smiles():
    assert parse_reaction_smiles(
        "CC.O>>CCO", ReactionFormat.STANDARD
    ) == ReactionEquation(["CC", "O"], [], ["CCO"])
    assert parse_reaction_smiles(
        "CC.O.[Na+]~[Cl-]>>CCO", ReactionFormat.STANDARD_WITH_TILDE
    ) == ReactionEquation(["CC", "O", "[Na+].[Cl-]"], [], ["CCO"])
    assert parse_reaction_smiles(
        "CC.O.[Na+].[Cl-]>>CCO |f:2.3|", ReactionFormat.EXTENDED
    ) == ReactionEquation(["CC", "O", "[Na+].[Cl-]"], [], ["CCO"])

    with pytest.raises(ValueError):
        _ = parse_reaction_smiles("CC.O>>CCO", "standard")  # type: ignore


def test_to_reaction_smiles():
    reaction_equation = ReactionEquation(["CC", "O", "[Na+].[Cl-]"], [], ["CCO"])
    assert (
        to_reaction_smiles(reaction_equation, ReactionFormat.STANDARD)
        == "CC.O.[Na+].[Cl-]>>CCO"
    )
    assert (
        to_reaction_smiles(reaction_equation, ReactionFormat.STANDARD_WITH_TILDE)
        == "CC.O.[Na+]~[Cl-]>>CCO"
    )
    assert (
        to_reaction_smiles(reaction_equation, ReactionFormat.EXTENDED)
        == "CC.O.[Na+].[Cl-]>>CCO |f:2.3|"
    )

    with pytest.raises(ValueError):
        _ = to_reaction_smiles(reaction_equation, "standard")  # type: ignore
