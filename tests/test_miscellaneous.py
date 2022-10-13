from collections import Counter

import pytest

from rxn.chemutils.exceptions import InvalidReactionSmiles, InvalidSmiles
from rxn.chemutils.miscellaneous import (
    apply_to_any_smiles,
    atom_type_counter,
    canonicalize_any,
    equivalent_smiles,
    is_valid_smiles,
    remove_chiral_centers,
    remove_double_bond_stereochemistry,
    sort_any,
)


def test_equivalent_smiles() -> None:
    assert equivalent_smiles("CCO", "C(C)O")
    assert not equivalent_smiles("CCO", "O(C)C")

    assert equivalent_smiles("CCO", "C(C)O", "OCC")
    assert not equivalent_smiles("CCO", "C(C)O", "O(C)C")
    assert not equivalent_smiles("CF=C", "F(C)C")

    # Should be false for invalid SMILES
    assert not equivalent_smiles("CCO", "", "O(C)C")

    # In case of invalid valence: not equivalent if check_valence flag is specified
    assert equivalent_smiles("CFC", "F(C)C")
    assert not equivalent_smiles("CFC", "F(C)C", check_valence=True)


def test_is_valid_smiles() -> None:
    assert is_valid_smiles("CCO")
    assert is_valid_smiles("C")
    assert is_valid_smiles("NOOOOO")
    assert is_valid_smiles("CC(CCC)")

    assert not is_valid_smiles("CFC")
    assert not is_valid_smiles("YES")
    assert not is_valid_smiles("Noo")
    assert not is_valid_smiles("CC(CC")

    assert is_valid_smiles("CFC", check_valence=False)
    assert not is_valid_smiles("YES", check_valence=False)
    assert not is_valid_smiles("CC(CC", check_valence=False)


def test_atom_type_counter() -> None:
    assert atom_type_counter("CCO") == Counter({"C": 2, "H": 6, "O": 1})
    assert atom_type_counter("CC=O") == Counter({"C": 2, "H": 4, "O": 1})
    assert atom_type_counter("[Na+].[Cl-]") == Counter({"Na": 1, "Cl": 1})

    # with radicals
    assert atom_type_counter("[C]CO") == Counter({"C": 2, "H": 3, "O": 1})

    with pytest.raises(InvalidSmiles):
        _ = atom_type_counter("CC(")


def test_remove_chiral_centers() -> None:
    input_expected_dict = {
        "O[C@](Br)(C)N": "O[C](Br)(C)N",
        "C[C@@](Br)(O)N": "C[C](Br)(O)N",
        "N[C@H](O)C": "N[CH](O)C",
        "N1[C@H](Cl)[C@@H](Cl)C(Cl)CC1": "N1[CH](Cl)[CH](Cl)C(Cl)CC1",
        "F/C=C/F": "F/C=C/F",
        "[N@+]": "[N+]",
        "[N@@+]": "[N+]",
        "[Si@@]": "[Si]",
        "[Si@H]": "[SiH]",
    }

    assert all(
        remove_chiral_centers(smi) == expected
        for smi, expected in input_expected_dict.items()
    )


def test_remove_double_bond_stereochemistry() -> None:
    input_expected_dict = {
        "C/C=C/C": "CC=CC",
        r"C/C=C\C": "CC=CC",
        r"C\C=C\C": "CC=CC",
        "F/C=C/C=C/C": "FC=CC=CC",
        "C1=CC=CC=C1": "C1=CC=CC=C1",
        "c1ccccc1": "c1ccccc1",
    }

    assert all(
        remove_double_bond_stereochemistry(smi) == expected
        for smi, expected in input_expected_dict.items()
    )


def test_apply_to_any_smiles() -> None:
    def dummy(smiles: str) -> str:
        return smiles + "0"

    # Single-component SMILES
    # Note: the first one is not considered to be a multi-component SMILES!
    assert apply_to_any_smiles("A.C.C.B", dummy) == "A.C.C.B0"
    assert apply_to_any_smiles("CBA", dummy) == "CBA0"

    # Multi-component SMILES (note that here, the fragments are not reordered)
    assert apply_to_any_smiles("A.D~C.B", dummy) == "A0.D~C0.B0"

    # reaction SMILES
    assert apply_to_any_smiles("B.A.E~D.A>>C.B", dummy) == "B0.A0.E~D0.A0>>C0.B0"
    assert apply_to_any_smiles("A.E.D.A>>C |f:1.2|", dummy) == "A0.A0.E.D0>>C0 |f:2.3|"


def test_canonicalize_any_on_molecule_smiles() -> None:
    assert canonicalize_any("CC(C)") == "CCC"
    assert canonicalize_any("[Na+].[Cl-]") == "[Cl-].[Na+]"
    assert canonicalize_any("CF(C)", check_valence=False) == "CFC"
    with pytest.raises(InvalidSmiles):
        _ = canonicalize_any("CF(C)")


def test_canonicalize_any_on_multicomponent_smiles() -> None:
    # Basic example
    assert canonicalize_any("CC(C).C(O)") == "CCC.CO"

    # Note: this will be assumed to be a single-compound SMILES -> reordered!
    assert canonicalize_any("C(O).CC(C)") == "CCC.CO"

    # Note: within the same compound: reordering due to compound canonicalization
    # Although "CCC" comes first alphabetically, it stays at the end.
    assert canonicalize_any("CO.C~O.CCC") == "CO.C~O.CCC"
    assert canonicalize_any("CO.O~C.CCC") == "CO.C~O.CCC"

    # Valence checking
    assert canonicalize_any("C(O)~CF(C)", check_valence=False) == "CFC~CO"
    with pytest.raises(InvalidSmiles):
        _ = canonicalize_any("CO.CFC")


def test_canonicalize_any_on_reaction_smiles() -> None:

    # Basic examples
    assert canonicalize_any("CC(C)>>C(O)") == "CCC>>CO"
    assert canonicalize_any("CC(C)>C(O)>C(O)") == "CCC>CO>CO"

    # Basic involving fragments
    assert canonicalize_any("CO.O~C>>C(O)") == "CO.C~O>>CO"
    assert canonicalize_any("CO.O.C>>C(O) |f:1.2|") == "CO.C.O>>CO |f:1.2|"
    # Note here the slight reordering - due to fragment parsing
    assert canonicalize_any("CC.CC.C.C>>CC |f:0.2|") == "CC.C.C.CC>>CC |f:2.3|"

    # Valence checking, and other potential exceptions
    assert canonicalize_any("CF(C)>>C(O)", check_valence=False) == "CFC>>CO"
    with pytest.raises(InvalidSmiles):
        _ = canonicalize_any("CFC>>CO")
    with pytest.raises(InvalidReactionSmiles):
        _ = canonicalize_any("CC>CO")
    with pytest.raises(InvalidReactionSmiles):
        _ = canonicalize_any("CC>>CC>>C(O)")


def test_sort_any() -> None:
    # Single-component SMILES
    assert sort_any("A.C.C.B") == "A.B.C.C"
    assert sort_any("CBA") == "CBA"

    # Multi-component SMILES (note that here, the fragments are not reordered)
    assert sort_any("A.D~C.B") == "A.B.D~C"

    # reaction SMILES
    assert sort_any("B.A.E~D.A>>C.B") == "A.A.B.E~D>>B.C"
    assert sort_any("B.A.E.D.A>>C.B |f:2.3|") == "A.A.B.E.D>>B.C |f:3.4|"
