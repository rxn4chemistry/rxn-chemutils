import pytest

from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.smiles_standardization import (
    standardize_molecules,
    standardize_smiles,
)


def test_standardize_smiles() -> None:
    smiles = "C(O)C"
    # case 1: canonicalization and sanitization
    assert standardize_smiles(smiles) == "CCO"
    # case 2: disabled canonicalization
    assert standardize_smiles(smiles, canonicalize=False) == "C(O)C"
    # case 3: disabled canonicalization with inchification (inherent canonicalization)
    assert standardize_smiles(smiles, canonicalize=False, inchify=True) == "CCO"
    # case 4: canonicalization different from inchification (tautomers interconversion)
    smiles = "CNC(=O)C"
    assert standardize_smiles(smiles) == "CNC(C)=O"
    assert standardize_smiles(smiles, inchify=True) == "CN=C(C)O"
    # case 5: inchification with metal disconnection
    smiles = "CCCC[Li]"
    assert standardize_smiles(smiles, inchify=True) == "[CH2]CCC.[Li]"
    # case 6: testing an invalid SMILES and error handling
    invalid_smiles = "C%5%%5"
    with pytest.raises(InvalidSmiles):
        standardize_smiles(invalid_smiles)


def test_standardize_molecules() -> None:
    # successful cases with different standardization flavours
    # case 1: default fragment bond
    molecules = "C(O)C.CCO.CC~C"
    assert standardize_molecules(molecules) == "CCO.CCO.C~CC"
    # case 2: custom fragment bond
    molecules = "C(O)C.CCO.CC|C"
    assert standardize_molecules(molecules, fragment_bond="|") == "CCO.CCO.C|CC"
    # case 3: molecule token delimiter
    molecules = "C(O)C.CCO.CC~C._C_"
    assert (
        standardize_molecules(molecules, molecule_token_delimiter="_")
        == "C.CCO.CCO.C~CC"
    )
    # case 4: molecule token delimiter with disabled ordering (canonicalization order the fragments)
    molecules = "C(O)C.CCO.CC~C._C_"
    assert (
        standardize_molecules(
            molecules, molecule_token_delimiter="_", ordered_precursors=False
        )
        == "CCO.CCO.C~CC.C"
    )
    # case 5: molecule token delimiter with disabled ordering and canonicalization
    molecules = "C(O)C.CCO.CC~C._C_"
    assert (
        standardize_molecules(
            molecules,
            canonicalize=False,
            molecule_token_delimiter="_",
            ordered_precursors=False,
        )
        == "C(O)C.CCO.CC~C.C"
    )
    # case 6: molecule token delimiter with disabled ordering and canonicalization, but enabled inchification
    molecules = "C(O)C.CCO.CC~C._C_"
    assert (
        standardize_molecules(
            molecules,
            canonicalize=False,
            inchify=True,
            molecule_token_delimiter="_",
            ordered_precursors=False,
        )
        == "CCO.CCO.C~CC.C"
    )
    # expected failures due to mismatch between the molecules string and the standard
    # case 7: unexpected fragment bond
    molecules = "C(O)C.CCO.CC|C"
    with pytest.raises(InvalidSmiles):
        standardize_molecules(molecules)
    # case 8: unexpected molecule token delimiter
    molecules = "C(O)C.CCO.CC~C._C_"
    with pytest.raises(InvalidSmiles):
        standardize_molecules(molecules)
