from pathlib import Path

import pytest

from rxn.chemutils.rdf import RdfParser, RdfReaction, ReactionSmilesExtractor

sample_rdf = Path(__file__).parent / "sample_with_unknown_structure.rdf"


@pytest.fixture
def reaction() -> RdfReaction:
    reactions = list(RdfParser(sample_rdf))
    return reactions[0]


def test_smiles_extraction_simple(reaction: RdfReaction) -> None:
    extractor = ReactionSmilesExtractor()
    # The unknown structure leads to a "*"
    assert (
        extractor.to_reaction_smiles(reaction)
        == "O=C1C=CC(=O)O1>CCOC(C)=O.*>O=C1OC(=O)C2C1C1C(=O)OC(=O)C21"
    )
