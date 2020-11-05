import fileinput

from rxn_chemutils.tokenization import tokenize_smiles


def tokenize() -> None:
    """
    Tokenize SMILES strings (molecules or reactions).

    The script will read strings either from stdin, or from a file given as the
    first argument (behavior from the fileinput package)."""

    for line in fileinput.input():
        smiles = line.strip()
        print(tokenize_smiles(smiles))


if __name__ == '__main__':
    tokenize()