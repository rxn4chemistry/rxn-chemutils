from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolToSmiles


class InvalidSmiles(ValueError):

    def __init__(self, smiles: str):
        super().__init__(f'"{smiles}" is not a valid SMILES string')


class SanitizationError(ValueError):

    def __init__(self, mol: Mol):
        message = 'Error when sanitizing RDKit Mol'
        try:
            smiles = MolToSmiles(mol)
            message += ': ' + smiles
        except Exception:
            pass
        super().__init__(message)
