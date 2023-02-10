import logging

import click
from rxn.utilities.files import (
    PathLike,
    dump_list_to_file,
    iterate_lines_from_file,
    raise_if_paths_are_identical,
)
from rxn.utilities.logging import setup_console_logger

from rxn.chemutils.multicomponent_smiles import (
    canonicalize_multicomponent_smiles,
    sort_multicomponent_smiles,
)
from rxn.chemutils.reaction_equation import canonicalize_compounds, sort_compounds
from rxn.chemutils.reaction_smiles import (
    determine_format,
    parse_reaction_smiles,
    to_reaction_smiles,
)

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def canonicalize_line(
    smiles_line: str, invalid_placeholder: str, sort_molecules: bool
) -> str:
    try:
        if ">" in smiles_line:
            # we have a reaction SMILES
            reaction_format = determine_format(smiles_line)
            reaction = parse_reaction_smiles(smiles_line, reaction_format)
            reaction = canonicalize_compounds(reaction)
            if sort_molecules:
                reaction = sort_compounds(reaction)
            return to_reaction_smiles(reaction, reaction_format)
        else:
            smiles = canonicalize_multicomponent_smiles(smiles_line, fragment_bond="~")
            if sort_molecules:
                smiles = sort_multicomponent_smiles(smiles)
            return smiles
    except Exception as e:
        logger.debug(f'Cannot canonicalize line "{smiles_line}": {e}')
        return invalid_placeholder


def canonicalize_file(
    input_file: PathLike,
    output_file: PathLike,
    invalid_placeholder: str = "",
    sort_molecules: bool = False,
) -> None:
    raise_if_paths_are_identical(input_file, output_file)
    logger.info(f'Canonicalizing file "{input_file}" -> "{output_file}".')

    # We formulate it as a generator, so that the file below is written directly
    canonical = (
        canonicalize_line(line, invalid_placeholder, sort_molecules)
        for line in iterate_lines_from_file(input_file)
    )

    dump_list_to_file(canonical, output_file)


@click.command()
@click.option("--input_file", "-i", required=True, help="File to canonicalize")
@click.option("--output_file", "-o", required=True, help="Canonicalized file")
@click.option(
    "--invalid_placeholder",
    default="",
    type=str,
    help="What to output when the canonicalization fails",
)
def main(input_file: str, output_file: str, invalid_placeholder: str) -> None:
    setup_console_logger()

    # Note: we put the actual code in a separate function, so that it can be
    # called also as a Python function.
    canonicalize_file(
        input_file=input_file,
        output_file=output_file,
        invalid_placeholder=invalid_placeholder,
    )


if __name__ == "__main__":
    main()
