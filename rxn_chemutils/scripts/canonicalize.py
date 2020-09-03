import click

from rxn_chemutils.utils import canonicalize_smiles, canonicalize_reaction_smiles


@click.command()
@click.option(
    '--input_file',
    '-i',
    required=True,
    help='File containing the SMILES strings to canonicalize (one per line)'
)
@click.option(
    '--output_file', '-o', required=True, help='File where to save the canonical SMILES strings'
)
def canonicalize(input_file: str, output_file: str) -> None:
    """Canonicalize SMILES strings (molecules or reactions) contained in a file"""

    with open(input_file, 'rt') as fin, open(output_file, 'wt') as fout:
        for line in fin:
            smiles = line.strip()

            if '>' in smiles:
                canonical_smiles = canonicalize_reaction_smiles(smiles)
            else:
                canonical_smiles = canonicalize_smiles(smiles)

            fout.write(f'{canonical_smiles}\n')


if __name__ == '__main__':
    canonicalize()
