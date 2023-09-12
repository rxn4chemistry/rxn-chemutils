import logging
from pathlib import Path
from typing import Iterable, Iterator

import click
from rxn.utilities.files import dump_list_to_file
from rxn.utilities.logging import setup_console_logger

from ..rdf import RdfReaction, ReactionSmilesExtractor, iterate_reactions_from_file

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@click.command()
@click.option(
    "-r",
    "--rdf_file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="RDF file (to read from).",
)
@click.option(
    "-s",
    "--smiles_file",
    type=click.Path(writable=True, path_type=Path),
    required=True,
    help="Where to save the reactions (in SMILES format, one per line).",
)
@click.option(
    "--fragment_bond",
    type=str,
    default="~",
    help="Fragment bond.",
)
def main(rdf_file: Path, smiles_file: Path, fragment_bond: str) -> None:
    """Convert a file of RDF reactions to SMILES format."""
    setup_console_logger()

    extractor = ReactionSmilesExtractor(fragment_bond=fragment_bond, sanitize=False)
    total_reactions = 0
    successful_reactions = 0

    def convert(rdfs: Iterable[RdfReaction]) -> Iterator[str]:
        """Convert reactions, ignoring the ones causing an error.

        Note: should be refactored to a class if the functionality is
        to be made a bit cleaner (avoid ``nonlocal``)."""
        nonlocal total_reactions
        nonlocal successful_reactions

        for rdf in rdfs:
            total_reactions += 1
            try:
                yield extractor.to_reaction_smiles(rdf)
                successful_reactions += 1
            except Exception as e:
                logger.warning(f"Cannot convert reaction: {e}")
                continue

    rdf_reactions = iterate_reactions_from_file(rdf_file)
    smiles_reactions = convert(rdf_reactions)
    dump_list_to_file(smiles_reactions, smiles_file)

    logger.info(
        f"Finished conversion. Successful: {successful_reactions} / {total_reactions}."
    )


if __name__ == "__main__":
    main()
