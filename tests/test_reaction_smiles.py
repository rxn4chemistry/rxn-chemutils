import pytest

from rxn.chemutils.reaction_equation import ReactionEquation
from rxn.chemutils.reaction_smiles import (
    ReactionFormat,
    canonicalize_reaction_smiles,
    determine_format,
    parse_any_reaction_smiles,
    parse_reaction_smiles,
    to_reaction_smiles,
)


def test_determine_format():
    assert determine_format("CC.O>>CCO") is ReactionFormat.STANDARD
    assert (
        determine_format("CC.O.[Na+]~[Cl-]>>CCO") is ReactionFormat.STANDARD_WITH_TILDE
    )
    assert determine_format("CC.O.[Na+].[Cl-]>>CCO |f:2.3|") is ReactionFormat.EXTENDED


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


def test_canonicalize_reaction_smiles():
    reaction_smiles = (
        "[CH:8]1=[CH:7][C:4](=[CH:3][C:2](=[CH:9]1)[Br:1])[CH:5]=[O:6].[CH2:10]([CH2:11]O)[OH:12]"
        ">CC1=CC=CC=C1.CC1=CC=C(C=C1)S(=O)(=O)O"
        ">[CH2:10]1[CH2:11][O:6][CH:5]([O:12]1)[C:4]2=[CH:3][C:2](=[CH:9][CH:8]=[CH:7]2)[Br:1]"
    )

    canonical = canonicalize_reaction_smiles(reaction_smiles)

    assert (
        canonical == "O=Cc1cccc(Br)c1.OCCO"
        ">Cc1ccccc1.Cc1ccc(S(=O)(=O)O)cc1"
        ">Brc1cccc(C2OCCO2)c1"
    )


def test_canonicalize_reaction_smiles_with_fragments():
    # No change to the fragment information should happen
    reaction_smiles = (
        "[CH3:1][C:2]([CH3:3])([CH3:4])[S@@:5](=[O:6])[NH2:7].[CH:14]1=[CH:13][C:11](=[CH:10][C:9](=[CH:15]1)[Br:8])[CH:12]=O"
        ">C(Cl)Cl.C(Cl)Cl.C(=O)([O-])[O-].[Cs+].[Cs+]"
        ">[CH3:1][C:2]([CH3:3])([CH3:4])[S@@:5](=[O:6])/[N:7]=[CH:12]/[C:11]1=[CH:10][C:9](=[CH:15][CH:14]=[CH:13]1)[Br:8]"
        " |f:4.5.6|"
    )

    canonical = canonicalize_reaction_smiles(reaction_smiles)

    assert (
        canonical == "CC(C)(C)[S@](N)=O.O=Cc1cccc(Br)c1"
        ">ClCCl.ClCCl.O=C([O-])[O-].[Cs+].[Cs+]"
        ">CC(C)(C)[S@@](=O)/N=C/c1cccc(Br)c1"
        " |f:4.5.6|"
    )
