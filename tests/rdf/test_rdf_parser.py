from pathlib import Path

from rxn.chemutils.rdf import RdfParser, ReactionSmilesExtractor

sample_rdf = Path(__file__).parent / "sample.rdf"

# structure present in the meta of the third reaction in the rdf
solvent_molstructure = """
  SOMEINFO          2D 1   1.00000     0.00000     0

  6  5  0  0  0  0  0  0  0  0999 V2000
  -11.0000   -3.3269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000   -3.3269    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000   -0.4425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    9.0000   -3.3269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    5.3269    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0000   -0.4425    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  6  1  0  0  0  0
  2  3  1  0  0  0  0
  2  6  1  0  0  0  0
  3  4  1  0  0  0  0
  3  5  2  0  0  0  0
M  END"""


def test_rdf_parser() -> None:
    rdf_parser = RdfParser(sample_rdf)
    reactions = list(rdf_parser.iter_reactions())

    # A few checks about the parsed reactions:

    # 1) number of parsed reactions
    assert len(reactions) == 3

    # 2) reaction number
    assert [reaction.reaction_index for reaction in reactions] == [1, 2, 6]

    # 3) ID2
    assert [reaction.meta["RXN:ID2"] for reaction in reactions] == [
        "RX0000100002000003",
        "RX0000100002000004",
        "RX0000100002000005",
    ]

    # 4) number of reactants and products for the first reaction
    assert len(reactions[0].reactants) == 2
    assert len(reactions[0].products) == 1

    # 5) embedded structure
    assert (
        reactions[2].meta["RXN:SOLVENT(1):MOL(1):MOLSTRUCTURE"] == solvent_molstructure
    )

    # 6) very last line of the file
    assert reactions[2].meta["RXN:PRODUCT(1):MOL_ID"] == "8"


def test_rdf_parser_extracted_smiles() -> None:
    rse = ReactionSmilesExtractor()
    assert [rse.to_reaction_smiles(reaction) for reaction in RdfParser(sample_rdf)] == [
        "O=C1C=CC(=O)N1.C#CCCCC>>CCCCC1=CC2C(=O)NC(=O)C12",
        "C=CCCCN1C(=O)C(C)=C(C)C1=O>>CC1=C(C)C(=O)N2CCCC2CC1=O",
        "O=C1C=CC(=O)O1>CCOC(C)=O>O=C1OC(=O)C2C1C1C(=O)OC(=O)C21",
    ]
