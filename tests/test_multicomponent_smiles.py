import pytest

from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.multicomponent_smiles import (
    apply_to_multicomponent_smiles,
    canonicalize_multicomponent_smiles,
    list_to_multicomponent_smiles,
    multicomponent_smiles_to_list,
    remove_duplicates_in_multicomponent_smiles,
    sort_multicomponent_smiles,
)


def test_multicomponent_smiles_to_list():
    assert multicomponent_smiles_to_list("") == []
    assert multicomponent_smiles_to_list("A.B.C") == ["A", "B", "C"]
    assert multicomponent_smiles_to_list("..A.B.C.") == ["A", "B", "C"]
    assert multicomponent_smiles_to_list("A~B.C") == ["A~B", "C"]
    assert multicomponent_smiles_to_list("A~B.C", fragment_bond="~") == ["A.B", "C"]
    assert multicomponent_smiles_to_list("A#B.C", fragment_bond="#") == ["A.B", "C"]


def test_list_to_multicomponent_smiles():
    assert list_to_multicomponent_smiles([]) == ""
    assert list_to_multicomponent_smiles(["A", "B", "C"]) == "A.B.C"
    assert list_to_multicomponent_smiles(["A.B", "C"]) == "A.B.C"
    assert list_to_multicomponent_smiles(["A.B", "C"], fragment_bond="~") == "A~B.C"
    assert list_to_multicomponent_smiles(["A.B", "C"], fragment_bond="#") == "A#B.C"


def test_apply_to_multicomponent_smiles():
    def dummy(smiles: str) -> str:
        """Dummy function that appends the number of fragments in a SMILES string to it."""
        return smiles + str(smiles.count(".") + 1)

    assert apply_to_multicomponent_smiles("A.B.CDE", dummy) == "A1.B1.CDE1"
    assert apply_to_multicomponent_smiles("A.B~C", dummy) == "A1.B~C1"
    assert (
        apply_to_multicomponent_smiles("A.B~C", dummy, fragment_bond="~") == "A1.B~C2"
    )
    assert (
        apply_to_multicomponent_smiles("ABCD.B~C~D", dummy, fragment_bond="~")
        == "ABCD1.B~C~D3"
    )


def test_canonicalize_multicomponent_smiles():
    # A few basic examples
    assert canonicalize_multicomponent_smiles("C(C)O") == "CCO"
    assert canonicalize_multicomponent_smiles("C.OCC.C(C)O") == "C.CCO.CCO"

    # Canonicalizing an empty string is valid - since this is a valid multi-component SMILES
    assert canonicalize_multicomponent_smiles("") == ""

    # with fragment bond; may fail if fragment bond not specified
    assert (
        canonicalize_multicomponent_smiles("C.[Na+]~[H-].O", fragment_bond="~")
        == "C.[H-]~[Na+].O"
    )
    with pytest.raises(InvalidSmiles):
        _ = canonicalize_multicomponent_smiles("C.[Na+]~[H-].O")

    # possibility to check valence or not
    assert (
        canonicalize_multicomponent_smiles("C.CF(C).O", check_valence=False)
        == "C.CFC.O"
    )
    with pytest.raises(InvalidSmiles):
        _ = canonicalize_multicomponent_smiles("C.CF(C).O")


def test_sort_multicomponent_smiles():
    assert sort_multicomponent_smiles("") == ""
    assert sort_multicomponent_smiles("B.D.C.A") == "A.B.C.D"
    assert sort_multicomponent_smiles("B.D.C~A") == "B.C~A.D"
    assert sort_multicomponent_smiles("NCC(N).CCO.COCO") == "CCO.COCO.NCC(N)"


def test_remove_duplicates_in_multicomponent_smiles():
    assert remove_duplicates_in_multicomponent_smiles("B.D.A.E") == "B.D.A.E"
    assert remove_duplicates_in_multicomponent_smiles("B.D.A.A.D.E") == "B.D.A.E"
    assert remove_duplicates_in_multicomponent_smiles("CC.C.O.N.O") == "CC.C.O.N"

    # Note: does not remove compounds that are the same when canonicalized!
    assert (
        remove_duplicates_in_multicomponent_smiles("CC.C(C).(C)(C)") == "CC.C(C).(C)(C)"
    )
