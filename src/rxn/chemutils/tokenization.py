import re
from typing import List

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


def tokenize_smiles(smiles: str) -> str:
    """
    Tokenize a SMILES molecule or reaction, and join the tokens with spaces.

    Args:
        smiles: SMILES string to tokenize, for instance 'CC(CO)=N>>CC(C=O)N'.

    Returns:
        SMILES string after tokenization, for instance 'C C ( C O ) = N >> C C ( C = O ) N'.
    """
    tokens = to_tokens(smiles)
    return " ".join(tokens)


def detokenize_smiles(tokenized_smiles: str) -> str:
    """
    Detokenize a tokenized SMILES string (that contains spaces between the characters).

    Args:
        tokenized_smiles: tokenized SMILES, for instance 'C C ( C O ) = N >> C C ( C = O ) N'

    Returns:
        SMILES after detokenization, for instance 'CC(CO)=N>>CC(C=O)N'
    """
    return tokenized_smiles.replace(" ", "")
