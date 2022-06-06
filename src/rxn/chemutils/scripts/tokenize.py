import sys
from typing import TextIO

import click

from rxn.chemutils.tokenization import tokenize_smiles


@click.command()
@click.argument("input_file", type=click.File(mode="r"), default=sys.stdin)
@click.argument("output_file", type=click.File(mode="w"), default=sys.stdout)
def main(input_file: TextIO, output_file: TextIO) -> None:
    """
    Tokenize SMILES strings (molecules or reactions).

    The script will read strings either from stdin, or from a file given as the
    first argument, and write to stdout, or from a file given as the second
    argument.
    """

    for line in input_file:
        smiles = line.strip()
        output_file.write(f"{tokenize_smiles(smiles)}\n")


if __name__ == "__main__":
    main()
