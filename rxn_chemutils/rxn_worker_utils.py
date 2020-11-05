"""Utilities used in the RXN forward/retrosynthesis workers."""
import sys
from rdkit import Chem
from io import StringIO
from loguru import logger
from typing import Optional

RXN_SMILES_SEPARATOR = '>>'


class RDKitError(ValueError):
    """Exception raised in RDKit."""

    def __init__(self, title: str, detail: str):
        """
        Initialize RDKitError.

        Args:
            title (str): title of the error.
            detail (str): decscription of the error.
        """
        self.type = 'RDKitError'
        self.title = title
        self.detail = detail


def standardize_smiles(
    smiles: str,
    canonicalize: bool = True,
    sanitize: bool = True,
    inchify: bool = False,
) -> str:
    """
    Ensure that a SMILES follows a desired standard.

    It allows canonicalization, sanitization and inchification keeping stereochemistry with isomericSmile=True.
    It can process multiple molecules separated by ".".
    Note that inchify set to True will also canonicalize the molecule.

    Args:
        smiles (str): SMILES representation of a molecule.
        canonicalize (bool): canonicalize SMILES. Defaults to True.
        sanitize (bool): sanitize SMILES. Defaults to True.
        inchify (bool): inchify the SMILES. Defaults to False.

    Raises:
        RDKitError: raised in case the SMILES does not represent a valid molecule
            or in the case the inchification failed if enabled.

    Returns:
        str: a SMILES following the desired standard.
    """
    molecule = Chem.MolFromSmiles(smiles, sanitize=sanitize)
    if molecule is None:
        raise RDKitError(
            'UnknownMolFromSmilesError', 'Invalid molecule for SMILES: {}'.format(smiles)
        )
    if inchify:
        inchi_string = Chem.MolToInchi(molecule)
        if inchi_string is None:
            logger.warning(
                'Inchification failure for SMILES: {}. Returning its canonical version.'
            )
            return Chem.MolToSmiles(molecule, isomericSmiles=True)
        else:
            # canonical set to True because we can't guarantee no canonicalization
            return Chem.MolToSmiles(Chem.MolFromInchi(inchi_string), canonical=True)
    if canonicalize:
        return Chem.MolToSmiles(molecule, isomericSmiles=True)
    else:
        return smiles


def standardize_molecules(
    molecules: str,
    canonicalize: bool = True,
    sanitize: bool = True,
    inchify: bool = False,
    fragment_bond: str = '~',
    ordered_precursors: bool = True,
    molecule_token_delimiter: Optional[str] = None
) -> str:
    """
    Ensure that a set of molecules represented by a string follows a desired standard.

    Args:
        molecules (str): molecules SMILES. Molecules can be separated via a ".".
            Fragments are supported with a custom `fragment_bond`.
        canonicalize (bool): canonicalize SMILES. Defaults to True.
        sanitize (bool): sanitize SMILES. Defaults to True.
        inchify (bool): inchify the SMILES. Defaults to False.
        fragment_bond (str): fragment bond. Defaults to '~'.
        ordered_precursors (bool): order precursors. Defaults to True.
        molecule_token_delimiter (Optional[str]): delimiter for big molecule tokens. Defaults to None


    Raises:
        RDKitError: errors in RDKit SMILES processing.
        Exception: generic uncaught error.

    Returns:
        str: standardized molecules.

    Examples:
        Standardize muliple molecules:

        >>> standardize_molecules('CCO.CC')
        'CC.CCO'

        Standardize multiple molecules including fragement information:

        >>> standardize_molecules('CCO.CC~C')
        'CCO.C~CC'
    """
    sio = sys.stderr = StringIO()
    Chem.WrapLogs()
    try:
        if molecule_token_delimiter is not None:
            molecules = molecules.replace(molecule_token_delimiter, '')
        if fragment_bond in molecules:
            standardized_molecules_list = [
                # make sure we remove the fragment to have valid SMILES
                standardize_smiles(
                    molecule.replace(fragment_bond, '.'),
                    canonicalize=canonicalize,
                    sanitize=sanitize,
                    inchify=inchify
                ).replace('.', fragment_bond) for molecule in molecules.split('.')
            ]
            if ordered_precursors:
                standardized_molecules_list = sorted(standardized_molecules_list)
            standardized_molecules = '.'.join(standardized_molecules_list)
        else:
            if ordered_precursors:
                # RDKit guarantees ordered precursors
                standardized_molecules = standardize_smiles(
                    molecules, canonicalize=canonicalize, sanitize=sanitize, inchify=inchify
                )
            else:
                standardized_molecules_list = [
                    standardize_smiles(
                        molecule, canonicalize=canonicalize, sanitize=sanitize, inchify=inchify
                    ) for molecule in molecules.split('.')
                ]
                standardized_molecules = '.'.join(standardized_molecules_list)
    except Exception:
        sio_str = str(sio.getvalue().strip())
        if 'Explicit valence for atom' in sio_str:
            raise RDKitError(
                'ExplicitValenceError', '{} - molecules: {}'.format(sio_str, molecules)
            )
        elif 'Could not sanitize molecule' in sio_str:
            raise RDKitError(
                'CouldNotSanitizeMoleculeError', '{} - molecules: {}'.format(sio_str, molecules)
            )
        elif 'unclosed ring' in sio_str:
            raise RDKitError('UnclosedRingError', '{} - molecules: {}'.format(sio_str, molecules))
        elif 'syntax error' in sio_str:
            raise RDKitError('SyntaxError', '{} - molecules: {}'.format(sio_str, molecules))
        elif "Can't kekulize mol" in '{} - molecules: {}'.format(sio_str, molecules):
            raise RDKitError('KekulizeMolError', sio_str)
        elif sio_str.startswith('RDKit ERROR:'):
            raise RDKitError(
                'UnknownStandardizationError', '{} - molecules: {}'.format(sio_str, molecules)
            )
        else:
            raise
    return standardized_molecules
