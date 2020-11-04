import pytest
from rxn_chemutils.rxn_worker_utils import RDKitError, standardize_smiles, standardize_molecules, tokenize


def test_rdkit_error():
    error = RDKitError('AnRDKitError', 'my message')
    assert error.type == 'RDKitError'
    assert error.title == 'AnRDKitError'
    assert error.detail == 'my message'
    with pytest.raises(ValueError):
        raise error


def test_standardize_smiles():
    smiles = 'C(O)C'
    assert standardize_smiles(smiles) == 'CCO'
    assert standardize_smiles(smiles, canonicalize=False) == 'C(O)C'
    assert standardize_smiles(smiles, canonicalize=False, inchify=True) == 'CCO'


def test_standardize_molecules():
    molecules = 'C(O)C.CCO.CC~C'
    assert standardize_molecules(molecules) == 'CCO.CCO.C~CC'
    molecules = 'C(O)C.CCO.CC|C'
    assert standardize_molecules(molecules, fragment_bond='|') == 'CCO.CCO.C|CC'
    molecules = 'C(O)C.CCO.CC~C._C_'
    assert standardize_molecules(molecules, molecule_token_delimiter='_') == 'C.CCO.CCO.C~CC'
    molecules = 'C(O)C.CCO.CC~C._C_'
    assert standardize_molecules(
        molecules, molecule_token_delimiter='_', ordered_precursors=False
    ) == 'CCO.CCO.C~CC.C'
    molecules = 'C(O)C.CCO.CC~C._C_'
    assert standardize_molecules(
        molecules, canonicalize=False, molecule_token_delimiter='_', ordered_precursors=False
    ) == 'C(O)C.CCO.CC~C.C'
    molecules = 'C(O)C.CCO.CC~C._C_'
    assert standardize_molecules(
        molecules,
        canonicalize=False,
        inchify=True,
        molecule_token_delimiter='_',
        ordered_precursors=False
    ) == 'CCO.CCO.C~CC.C'


def test_tokenize():
    smiles = 'COCO.OCC>>NCOC'
    expected = 'C O C O . O C C >> N C O C'.split()
    assert tokenize(smiles) == expected
