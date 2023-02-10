import click
from rxn.utilities.files import load_list_from_file

from rxn.chemutils.reaction_combiner import ReactionCombiner
from rxn.chemutils.reaction_smiles import ReactionFormat


@click.command()
@click.argument("fragments_1_file", type=click.Path(exists=True, dir_okay=False))
@click.argument("fragments_2_file", type=click.Path(exists=True, dir_okay=False))
@click.option("--standardize/--no-standardize", default=False)
@click.option(
    "--reaction_format",
    type=click.Choice(
        ["standard", "standard_with_tilde", "extended"], case_sensitive=False
    ),
    default="extended",
)
def main(
    fragments_1_file: str,
    fragments_2_file: str,
    standardize: bool,
    reaction_format: str,
) -> None:
    """Combine precursors and products (or two sets of partial reactions, and
    write the reactions into std output.

    If one of both file sizes is a multiple of the other, we consider it to
    have been generated via a "top-N" prediction.
    """
    fragments_1 = load_list_from_file(fragments_1_file)
    fragments_2 = load_list_from_file(fragments_2_file)

    combiner = ReactionCombiner(
        standardize=standardize,
        reaction_format=ReactionFormat.from_string(reaction_format),
    )

    for reaction_smiles in combiner.combine(fragments_1, fragments_2):
        print(reaction_smiles)


if __name__ == "__main__":
    main()
