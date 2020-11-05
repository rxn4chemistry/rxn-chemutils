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
    invalid_smiles = 'C%5%%5'
    with pytest.raises(RDKitError):
        try:
            standardize_smiles(invalid_smiles)
        except RDKitError as exception:
            assert exception.title == 'UnknownMolFromSmilesError'
            assert exception.detail == 'Invalid molecule for SMILES: {}'.format(invalid_smiles)
            raise


def test_standardize_molecules():
    # successful cases with different standardization flavours
    # case 1: default fragment bond
    molecules = 'C(O)C.CCO.CC~C'
    assert standardize_molecules(molecules) == 'CCO.CCO.C~CC'
    # case 2: custom fragment bond
    molecules = 'C(O)C.CCO.CC|C'
    assert standardize_molecules(molecules, fragment_bond='|') == 'CCO.CCO.C|CC'
    # case 3: molecule token delimiter
    molecules = 'C(O)C.CCO.CC~C._C_'
    assert standardize_molecules(molecules, molecule_token_delimiter='_') == 'C.CCO.CCO.C~CC'
    # case 4: molecule token delimiter with disabled ordering (canonicalization order the fragments)
    molecules = 'C(O)C.CCO.CC~C._C_'
    assert standardize_molecules(
        molecules, molecule_token_delimiter='_', ordered_precursors=False
    ) == 'CCO.CCO.C~CC.C'
    # case 5: molecule token delimiter with disabled ordering and canonicalization
    molecules = 'C(O)C.CCO.CC~C._C_'
    assert standardize_molecules(
        molecules, canonicalize=False, molecule_token_delimiter='_', ordered_precursors=False
    ) == 'C(O)C.CCO.CC~C.C'
    # case 6: molecule token delimiter with disabled ordering and canonicalization, but enabled inchification
    molecules = 'C(O)C.CCO.CC~C._C_'
    assert standardize_molecules(
        molecules,
        canonicalize=False,
        inchify=True,
        molecule_token_delimiter='_',
        ordered_precursors=False
    ) == 'CCO.CCO.C~CC.C'
    # expected failures due to mismatch between the molecules string and the standard
    # case 7: unexpected fragment bond
    molecules = 'C(O)C.CCO.CC|C'
    with pytest.raises(RDKitError):
        try:
            standardize_molecules(molecules)
        except RDKitError as exception:
            assert exception.title == 'UnknownMolFromSmilesError'
            assert exception.detail == 'Invalid molecule for SMILES: {}'.format(molecules)
            raise
    # case 8: unexpected molecule token delimiter
    molecules = 'C(O)C.CCO.CC~C._C_'
    with pytest.raises(RDKitError):
        try:
            standardize_molecules(molecules)
        except RDKitError as exception:
            assert exception.title == 'UnknownMolFromSmilesError'
            # since we have a fragment bond we expect the single molecules to be processed
            # one by one, hence the error should be raised only in the one with the unstripped
            # molecule token delimiters.
            assert exception.detail == 'Invalid molecule for SMILES: _C_'
            raise


def test_tokenize():
    smiles = 'COCO.OCC>>NCOC'
    expected = 'C O C O . O C C >> N C O C'.split()
    assert tokenize(smiles) == expected
