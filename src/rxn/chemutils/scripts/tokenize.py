import sys
from typing import Optional, TextIO

import click
from rxn.utilities.logging import setup_console_logger

from rxn.chemutils.tokenization import tokenize_smiles


@click.command()
@click.argument("input_file", type=click.File(mode="r"), default=sys.stdin)
@click.argument("output_file", type=click.File(mode="w"), default=sys.stdout)
@click.option(
    "--fallback_value",
    type=str,
    help="Placeholder for strings that cannot be tokenized. By default, an exception is raised.",
)
def main(
    input_file: TextIO, output_file: TextIO, fallback_value: Optional[str]
) -> None:
    """
    Tokenize SMILES strings (molecules or reactions).

    The script will read strings either from stdin, or from a file given as the
    first argument, and write to stdout, or from a file given as the second
    argument.
    """
    setup_console_logger()

    for line in input_file:
        smiles = line.strip()
        output_file.write(f"{tokenize_smiles(smiles, fallback_value=fallback_value)}\n")


if __name__ == "__main__":
    main()
