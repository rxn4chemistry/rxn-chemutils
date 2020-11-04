from typing import Sequence

from rdkit.Chem.rdchem import Mol


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


def remove_atom_mapping(mols: Sequence[Mol]) -> None:
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                atom.ClearProp('molAtomMapNumber')
