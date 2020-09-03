import re
from typing import List


def to_tokens(smiles: str) -> List[str]:
    """
    Tokenize a SMILES molecule or reaction into a list of tokens.

    Args:
        smiles: SMILES string to tokenize.

    Returns:
        List of tokens (give back the original SMILES string if appended).
    """
    pattern = r"(\%\([0-9]{3}\)|\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\|\/|:|~|@|\?|>>?|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smiles)]

    assert smiles == ''.join(tokens), f'SMI: {smiles} =!= JOIN: {"".join(tokens)}'

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
    return ' '.join(tokens)


def detokenize_smiles(tokenized_smiles: str) -> str:
    """
    Detokenize a tokenized SMILES string (that contains spaces between the characters).

    Args:
        tokenized_smiles: tokenized SMILES, for instance 'C C ( C O ) = N >> C C ( C = O ) N'

    Returns:
        SMILES after detokenization, for instance 'CC(CO)=N>>CC(C=O)N'
    """
    return tokenized_smiles.replace(' ', '')
