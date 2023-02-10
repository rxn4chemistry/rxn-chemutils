import pytest
from rxn.utilities.files import (
    dump_list_to_file,
    load_list_from_file,
    named_temporary_path,
)

from rxn.chemutils.exceptions import UnclearWhetherTokenized
from rxn.chemutils.tokenization import (
    TokenizationError,
    detokenize_file,
    detokenize_smiles,
    ensure_tokenized_file,
    file_is_tokenized,
    string_is_tokenized,
    to_tokens,
    tokenize_file,
    tokenize_smiles,
)


def test_tokenization_error() -> None:
    error = TokenizationError("ATokenizationError", "my message")
    assert error.type == "TokenizationError"
    assert error.title == "ATokenizationError"
    assert error.detail == "my message"
    with pytest.raises(ValueError):
        raise error


def test_to_tokens() -> None:
    smiles = r"C\C=C\C"
    groundtruth = ["C", "\\", "C", "=", "C", "\\", "C"]
    tokens = to_tokens(smiles)
    assert tokens == groundtruth
    smiles = "CCZCC"
    with pytest.raises(TokenizationError):
        try:
            to_tokens(smiles)
        except TokenizationError as exception:
            assert exception.title == "SmilesJoinedTokensMismatch"
            assert exception.detail == 'SMILES="CCZCC" != joined_tokens="CCCC"'
            raise


def test_tokenize_simple_smiles() -> None:
    # For this simple SMILES, 1 letter per token
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    expected = "C N 1 C = N C 2 = C 1 C ( = O ) N ( C ( = O ) N 2 C ) C"
    assert tokenize_smiles(smiles) == expected


def test_tokenize_complex_smiles() -> None:
    # Some tokens have more than one letter
    smiles = "C([N+](=O)[O-])1=CC(Br)=CC=C1"
    expected = "C ( [N+] ( = O ) [O-] ) 1 = C C ( Br ) = C C = C 1"
    assert tokenize_smiles(smiles) == expected


def test_tokenize_two_digit_cycle() -> None:
    smiles = "C%42CCCCC%42"
    expected = "C %42 C C C C C %42"
    assert tokenize_smiles(smiles) == expected


def test_tokenize_three_digit_cycle() -> None:
    smiles = "C%(432)CCCCC%(432)"
    expected = "C %(432) C C C C C %(432)"
    assert tokenize_smiles(smiles) == expected


def test_tokenize_reaction_smiles() -> None:
    smiles = "COCO.OCC>>NCOC"
    expected = "C O C O . O C C >> N C O C"
    assert tokenize_smiles(smiles) == expected


def test_tokenize_reaction_smiles_with_fragment_bond() -> None:
    smiles = "COCO.[Na+]~[OH-].OCC>O.C>NCOC"
    expected = "C O C O . [Na+] ~ [OH-] . O C C > O . C > N C O C"
    assert tokenize_smiles(smiles) == expected


def test_detokenize() -> None:
    # Basically, detokenizing simply removes the spaces
    tokenized = "some random string also not S M I L E S ."
    expected = "somerandomstringalsonotSMILES."
    assert detokenize_smiles(tokenized) == expected


def test_tokenization_roundtrip() -> None:
    smiles_list = [
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # simple SMILES
        "C([N+](=O)[O-])1=CC(Br)=CC=C1",  # more complex SMILES
        "C%42CCCCC%42",  # SMILES with two-digit cycle
        "C%(432)CCCCC%(432)",  # SMILES with three-digit cycle
        "COCO.OCC>>NCOC",  # reaction SMILES
        "COCO.[Na+]~[OH-].OCC>O.C>NCOC",  # reaction SMILES with fragment bonds
    ]

    for smiles in smiles_list:
        assert detokenize_smiles(tokenize_smiles(smiles)) == smiles


def test_string_is_tokenized():
    assert string_is_tokenized("C C O . [Na+]")
    assert string_is_tokenized("C C O . [Na+] >> C ( O ) Cl")
    assert not string_is_tokenized("C C O . [Na+] >> C (O) Cl")
    assert not string_is_tokenized("CCO")

    # Special cases - empty strings and strings with single tokens
    for string in ["", "C", ">>", "[Na]"]:
        with pytest.raises(UnclearWhetherTokenized):
            _ = string_is_tokenized(string)

    # Tokenization errors are being propagated
    with pytest.raises(TokenizationError):
        _ = string_is_tokenized("I N V A L I D")


def test_file_is_tokenized():
    # Basic tokenized example
    with named_temporary_path() as path:
        dump_list_to_file(["C C O >> C C O", "C C . C"], path)
        assert file_is_tokenized(path)

    # Basic non-tokenized example
    with named_temporary_path() as path:
        dump_list_to_file(["CCO>>CCO", "CC.C"], path)
        assert not file_is_tokenized(path)

    # Only checks the first line - returns True even if the second line is not tokenized
    with named_temporary_path() as path:
        dump_list_to_file(["C C O >> C C O", "CC.C"], path)
        assert file_is_tokenized(path)

    # empty file
    with named_temporary_path() as path:
        dump_list_to_file([], path)
        with pytest.raises(RuntimeError):
            _ = file_is_tokenized(path)

    # Invalid SMILES
    with named_temporary_path() as path:
        dump_list_to_file(["I N V A L I D", "CC.C"], path)
        with pytest.raises(TokenizationError):
            _ = file_is_tokenized(path)

    # Empty first line - needs to check the second line!
    with named_temporary_path() as path:
        dump_list_to_file(["", "C C O >> C C O"], path)
        assert file_is_tokenized(path)
    with named_temporary_path() as path:
        dump_list_to_file(["", "CCO>>CCO"], path)
        assert not file_is_tokenized(path)

    # First line has one single token - needs to check the second line!
    with named_temporary_path() as path:
        dump_list_to_file([">>", "C C O >> C C O"], path)
        assert file_is_tokenized(path)
    with named_temporary_path() as path:
        dump_list_to_file([">>", "CCO>>CCO"], path)
        assert not file_is_tokenized(path)


def test_tokenize_file():
    with named_temporary_path() as f_in, named_temporary_path() as f_out:
        # Original content
        original = ["CCO>>CCO", "CC.C", "INVALID", "C(NCC)[S]OC"]
        dump_list_to_file(original, f_in)

        # Expected (tokenized) content
        placeholder = "ERROR"
        tokenized = ["C C O >> C C O", "C C . C", placeholder, "C ( N C C ) [S] O C"]

        tokenize_file(f_in, f_out, invalid_placeholder=placeholder)

        assert load_list_from_file(f_out) == tokenized


def test_detokenize_file():
    with named_temporary_path() as f_in, named_temporary_path() as f_out:
        # Original (tokenized) content
        original = ["C C O >> C C O", "C C . C", "C ( N C C ) [S] O C"]
        dump_list_to_file(original, f_in)

        # Expected (detokenized) content
        detokenized = ["CCO>>CCO", "CC.C", "C(NCC)[S]OC"]

        detokenize_file(f_in, f_out)

        assert load_list_from_file(f_out) == detokenized


def test_ensure_tokenized_file():
    with named_temporary_path() as temporary_path:
        temporary_path.mkdir()

        # prepare filenames
        postfix = ".tknz"
        already_tokenized_file = str(temporary_path / "a.txt")
        not_tokenized_file = str(temporary_path / "b.txt")
        updated_tokenized_file = str(temporary_path / "b.txt") + postfix

        # contents (original and expected)
        placeholder = "error"
        tokenized = ["C C O >> C C O", "C C . C", "C ( N C C ) [S] O"]
        not_tokenized = ["CCO>>CCO", "CC.C", "INVALID", "C(N)[S]O"]
        after_tokenization = ["C C O >> C C O", "C C . C", placeholder, "C ( N ) [S] O"]

        # Put into files
        dump_list_to_file(tokenized, already_tokenized_file)
        dump_list_to_file(not_tokenized, not_tokenized_file)

        # ensure for already tokenized - the original unchanged file can be used
        result = ensure_tokenized_file(
            already_tokenized_file, postfix=postfix, invalid_placeholder=placeholder
        )
        assert result == already_tokenized_file
        assert load_list_from_file(result) == tokenized

        # ensure for non-tokenized - a new file was created with tokenized strings
        result = ensure_tokenized_file(
            not_tokenized_file, postfix=postfix, invalid_placeholder=placeholder
        )
        assert result == updated_tokenized_file
        assert load_list_from_file(updated_tokenized_file) == after_tokenization
