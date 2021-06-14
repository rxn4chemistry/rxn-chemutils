# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolToSmiles


class InvalidSmiles(ValueError):

    def __init__(self, smiles: str):
        super().__init__(f'"{smiles}" is not a valid SMILES string')


class InvalidMdl(ValueError):
    """
    Exception raised when converting invalid MDL Mol strings.
    """

    def __init__(self, mdl: str):
        super().__init__(f'The following MDL string cannot be converted: {mdl}')


class SanitizationError(ValueError):

    def __init__(self, mol: Mol):
        message = 'Error when sanitizing RDKit Mol'
        try:
            smiles = MolToSmiles(mol)
            message += ': ' + smiles
        except Exception:
            pass
        super().__init__(message)
