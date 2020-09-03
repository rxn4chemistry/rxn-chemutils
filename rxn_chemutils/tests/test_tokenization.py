from rxn_chemutils.tokenization import tokenize_smiles, detokenize_smiles


def test_tokenize_simple_smiles():
    # For this simple SMILES, 1 letter per token
    smiles = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
    expected = 'C N 1 C = N C 2 = C 1 C ( = O ) N ( C ( = O ) N 2 C ) C'
    assert tokenize_smiles(smiles) == expected


def test_tokenize_complex_smiles():
    # Some tokens have more than one letter
    smiles = 'C([N+](=O)[O-])1=CC(Br)=CC=C1'
    expected = 'C ( [N+] ( = O ) [O-] ) 1 = C C ( Br ) = C C = C 1'
    assert tokenize_smiles(smiles) == expected


def test_tokenize_two_digit_cycle():
    smiles = 'C%42CCCCC%42'
    expected = 'C %42 C C C C C %42'
    assert tokenize_smiles(smiles) == expected


def test_tokenize_three_digit_cycle():
    smiles = 'C%(432)CCCCC%(432)'
    expected = 'C %(432) C C C C C %(432)'
    assert tokenize_smiles(smiles) == expected


def test_tokenize_reaction_smiles():
    smiles = 'COCO.OCC>>NCOC'
    expected = 'C O C O . O C C >> N C O C'
    assert tokenize_smiles(smiles) == expected


def test_tokenize_reaction_smiles_with_fragment_bond():
    smiles = 'COCO.[Na+]~[OH-].OCC>O.C>NCOC'
    expected = 'C O C O . [Na+] ~ [OH-] . O C C > O . C > N C O C'
    assert tokenize_smiles(smiles) == expected


def test_detokenize():
    # Basically, detokenizing simply removes the spaces
    tokenized = 'some random string also not S M I L E S .'
    expected = 'somerandomstringalsonotSMILES.'
    assert detokenize_smiles(tokenized) == expected


def test_tokenization_roundtrip():
    smiles_list = [
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # simple SMILES
        'C([N+](=O)[O-])1=CC(Br)=CC=C1',  # more complex SMILES
        'C%42CCCCC%42',  # SMILES with two-digit cycle
        'C%(432)CCCCC%(432)',  # SMILES with three-digit cycle
        'COCO.OCC>>NCOC',  # reaction SMILES
        'COCO.[Na+]~[OH-].OCC>O.C>NCOC',  # reaction SMILES with fragment bonds
    ]

    for smiles in smiles_list:
        assert detokenize_smiles(tokenize_smiles(smiles)) == smiles
