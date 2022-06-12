from typing import Iterable

from rdkit.Chem import CombineMols
from rdkit.Chem.rdchem import Mol


def combine_mols(mols: Iterable[Mol]) -> Mol:
    """Combine multiple RDKit Mols into one.

    RDKit has a function for that but it only supports two mols.

    Args:
        mols: RDKit Mols to combine.
    """

    resulting_mol = Mol()
    for mol in mols:
        resulting_mol = CombineMols(resulting_mol, mol)
    return resulting_mol
