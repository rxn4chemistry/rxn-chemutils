import sys
from typing import Optional, TextIO

import click
from rxn.utilities.logging import setup_console_logger

from rxn.chemutils.miscellaneous import canonicalize_any


@click.command()
@click.argument("input_file", type=click.File(mode="r"), default=sys.stdin)
@click.argument("output_file", type=click.File(mode="w"), default=sys.stdout)
@click.option(
    "--invalid_placeholder",
    type=str,
    help=(
        "If specified, the given value will be the output for invalid input SMILES. "
        "By default, an exception is raised in such cases."
    ),
)
@click.option(
    "--sort_compounds",
    "-s",
    is_flag=True,
    help="If specified, the compounds will be sorted after canonicalization.",
)
def main(
    input_file: TextIO,
    output_file: TextIO,
    invalid_placeholder: Optional[str],
    sort_compounds: bool,
) -> None:
    """
    Canonicalize SMILES strings (molecules, sets of molecules, or reactions).

    The script will read strings either from stdin, or from a file given as the
    first argument, and write to stdout, or from a file given as the second
    argument.
    """
    setup_console_logger()

    for line in input_file:
        smiles = line.strip()

        # Canonicalize the SMILES, handle exception if needed
        try:
            canonical = canonicalize_any(smiles, sort_molecules=sort_compounds)
        except Exception:
            if invalid_placeholder is None:
                raise
            canonical = invalid_placeholder

        output_file.write(f"{canonical}\n")


if __name__ == "__main__":
    main()
