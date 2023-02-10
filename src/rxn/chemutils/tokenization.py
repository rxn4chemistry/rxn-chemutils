import logging
import re
import shutil
from typing import List, Optional

from rxn.utilities.files import (
    PathLike,
    dump_list_to_file,
    iterate_lines_from_file,
    raise_if_paths_are_identical,
)

from .exceptions import UnclearWhetherTokenized

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


SMILES_TOKENIZER_PATTERN = r"(\%\([0-9]{3}\)|\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\||\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"
SMILES_REGEX = re.compile(SMILES_TOKENIZER_PATTERN)


class TokenizationError(ValueError):
    """Exception raised in RDKit."""

    def __init__(self, title: str, detail: str):
        """
        Initialize TokenizationError.

        Args:
            title (str): title of the error.
            detail (str): decscription of the error.
        """
        self.type = "TokenizationError"
        self.title = title
        self.detail = detail


def to_tokens(smiles: str) -> List[str]:
    """
    Tokenize a SMILES molecule or reaction into a list of tokens.

    Args:
        smiles: SMILES string to tokenize.

    Raises:
        TokenizationError: in case of mismatch between the SMILES and the joined tokens.

    Returns:
        List of tokens (give back the original SMILES string if appended).
    """
    tokens = [token for token in SMILES_REGEX.findall(smiles)]

    if smiles != "".join(tokens):
        raise TokenizationError(
            "SmilesJoinedTokensMismatch",
            f'SMILES="{smiles}" != joined_tokens="{"".join(tokens)}"',
        )

    return tokens


def tokenize_smiles(smiles: str, fallback_value: Optional[str] = None) -> str:
    """
    Tokenize a SMILES molecule or reaction, and join the tokens with spaces.

    Args:
        smiles: SMILES string to tokenize, for instance 'CC(CO)=N>>CC(C=O)N'.
        fallback_value: what value to returns when the tokenization is unsuccessful.
            Default: no fallback, will propagate the TokenizationError exception.

    Returns:
        SMILES string after tokenization, for instance 'C C ( C O ) = N >> C C ( C = O ) N'.
    """
    try:
        tokens = to_tokens(smiles)
        return " ".join(tokens)
    except TokenizationError:
        if fallback_value is not None:
            logger.debug(f'Error when tokenizing "{smiles}"')
            return fallback_value
        raise


def detokenize_smiles(tokenized_smiles: str) -> str:
    """
    Detokenize a tokenized SMILES string (that contains spaces between the characters).

    Args:
        tokenized_smiles: tokenized SMILES, for instance 'C C ( C O ) = N >> C C ( C = O ) N'

    Returns:
        SMILES after detokenization, for instance 'CC(CO)=N>>CC(C=O)N'
    """
    return tokenized_smiles.replace(" ", "")


def string_is_tokenized(smiles_line: str) -> bool:
    """
    Whether a string is a tokenized SMILES or not.

    Args:
        smiles_line: string to inspect

    Raises:
        ValueError: if not possible to determine whether tokenized or not
        TokenizationError: propagated directly from tokenize_smiles()
    """
    detokenized = detokenize_smiles(smiles_line)
    tokens = to_tokens(detokenized)
    if len(tokens) < 2:
        raise UnclearWhetherTokenized(smiles_line)
    return " ".join(tokens) == smiles_line


def tokenize_file(
    input_file: PathLike, output_file: PathLike, fallback_value: str = ""
) -> None:
    """
    Tokenize a file containing SMILES strings.

    Args:
        input_file: file to tokenize.
        output_file: where to save the tokenized file.
        fallback_value: placeholder for strings that cannot be tokenized.
    """
    raise_if_paths_are_identical(input_file, output_file)
    logger.info(f'Tokenizing "{input_file}" -> "{output_file}".')

    tokenized = (
        tokenize_smiles(line, fallback_value)
        for line in iterate_lines_from_file(input_file)
    )

    dump_list_to_file(tokenized, output_file)


def detokenize_file(
    input_file: PathLike,
    output_file: PathLike,
) -> None:
    raise_if_paths_are_identical(input_file, output_file)
    logger.info(f'Detokenizing "{input_file}" -> "{output_file}".')

    detokenized = (
        detokenize_smiles(line) for line in iterate_lines_from_file(input_file)
    )
    dump_list_to_file(detokenized, output_file)


def ensure_tokenized_file(
    file: PathLike, postfix: str = ".tokenized", fallback_value: str = ""
) -> str:
    """
    Ensure that a file is tokenized: do nothing if the file is already tokenized, create
    a tokenized copy otherwise.

    Args:
        file: path to the file that we want to ensure is tokenized.
        postfix: postfix to add to the tokenized copy (if applicable).
        fallback_value: placeholder for strings that cannot be tokenized (if applicable).

    Returns:
        The path to the tokenized file (original path, or path to new file).
    """
    if file_is_tokenized(file):
        return str(file)

    tokenized_copy = str(file) + postfix
    tokenize_file(file, tokenized_copy, fallback_value=fallback_value)
    return tokenized_copy


def file_is_tokenized(filepath: PathLike) -> bool:
    """
    Whether a file contains tokenized SMILES or not.

    By default, this looks at the first non-empty line of the file only!

    Raises:
        TokenizationError: propagated from tokenize_smiles()
        RuntimeError: for empty files or files with empty lines only.

    Args:
        filepath: path to the file.
    """
    # Iterative formulation in case the first line(s) of the file don't make it
    # clear whether tokenized or not.
    for line in iterate_lines_from_file(filepath):
        try:
            return string_is_tokenized(line)
        except UnclearWhetherTokenized:
            continue
    raise RuntimeError(
        f'Could not determine whether "{filepath}" is tokenized: empty lines only.'
    )


def copy_as_detokenized(src: PathLike, dest: PathLike) -> None:
    """
    Copy a source file to a destination, while making sure that it is not tokenized.
    """
    if file_is_tokenized(src):
        logger.info(f'Copying and detokenizing "{src}" -> "{dest}".')
        detokenize_file(src, dest)
    else:
        logger.info(f'Copying "{src}" -> "{dest}".')
        shutil.copy(src, dest)
