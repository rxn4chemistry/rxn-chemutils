from typing import Optional

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolToSmiles


class InvalidSmiles(ValueError):
    def __init__(self, smiles: str, msg: Optional[str] = None):
        if msg is None:
            msg = f'"{smiles}" is not a valid SMILES string'
        super().__init__(msg)
        self.smiles = smiles


class InvalidReactionSmiles(InvalidSmiles):
    def __init__(self, reaction_smiles: str, msg: Optional[str] = None):
        if msg is None:
            msg = f'"{reaction_smiles}" is not a valid reaction SMILES string'
        super().__init__(smiles=reaction_smiles, msg=msg)


class InvalidMdl(ValueError):
    """
    Exception raised when converting invalid MDL Mol strings.
    """

    def __init__(self, mdl: str):
        super().__init__(f"The following MDL string cannot be converted: {mdl}")
        self.mdl = mdl


class SanitizationError(ValueError):
    def __init__(self, mol: Mol):
        message = "Error when sanitizing RDKit Mol"
        try:
            smiles = MolToSmiles(mol)
            message += ": " + smiles
        except Exception:
            pass
        super().__init__(message)
