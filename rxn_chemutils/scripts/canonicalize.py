import fileinput

from rxn_chemutils.conversion import canonicalize_smiles, canonicalize_reaction_smiles


def canonicalize() -> None:
    """
    Canonicalize SMILES strings (molecules or reactions).

    The script will read strings either from stdin, or from a file given as the
    first argument (behavior from the fileinput package)."""

    for line in fileinput.input():
        smiles = line.strip()

        if '>' in smiles:
            canonical_smiles = canonicalize_reaction_smiles(smiles)
        else:
            canonical_smiles = canonicalize_smiles(smiles)

        print(canonical_smiles)


if __name__ == '__main__':
    canonicalize()
