import click
from rxn.utilities.files import load_list_from_file

from rxn.chemutils.reaction_combiner import ReactionCombiner
from rxn.chemutils.reaction_smiles import ReactionFormat


@click.command()
@click.argument("precursors_file")
@click.argument("products_file")
@click.option("--no-standardize", is_flag=True)
@click.option(
    "--reaction_format",
    type=click.Choice(
        ["standard", "standard_with_tilde", "extended"], case_sensitive=False
    ),
    default="extended",
)
def main(
    precursors_file: str, products_file: str, no_standardize: bool, reaction_format: str
) -> None:
    """Combine precursors and products, and write the reactions into std output.

    If one of both file sizes is a multiple of the other, we consider it to
    have been generated via a "top-N" prediction.
    """
    precursors = load_list_from_file(precursors_file)
    products = load_list_from_file(products_file)

    combiner = ReactionCombiner(
        standardize=not no_standardize,
        reaction_format=ReactionFormat.from_string(reaction_format),
    )

    for reaction_smiles in combiner.combine(precursors, products):
        print(reaction_smiles)


if __name__ == "__main__":
    main()
