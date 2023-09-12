from rxn.chemutils.rdf import RdfReaction, ReactionSmilesExtractor

dummy_molstructure = """
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


def test_to_reaction_smiles() -> None:
    # In an earlier version, the reaction SMILES changed every time
    # to_reaction_smiles() was called.
    reaction = RdfReaction(
        [dummy_molstructure],
        [],
        [dummy_molstructure],
        meta={
            "RXN:SOME_KEY": "dummy",
            "RXN:CATALYST(6):MOL(1):MOLSTRUCTURE": dummy_molstructure,
        },
        reaction_index=22,
    )
    rse = ReactionSmilesExtractor()
    assert rse.to_reaction_smiles(reaction) == "CCOC(C)=O>CCOC(C)=O>CCOC(C)=O"
    # call it again
    assert rse.to_reaction_smiles(reaction) == "CCOC(C)=O>CCOC(C)=O>CCOC(C)=O"
