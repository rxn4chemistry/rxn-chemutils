from typing import Sequence

from rdkit.Chem.rdchem import Mol


def clear_atom_mapping(mols: Sequence[Mol]) -> None:
    """
    Remove the mapping information in RDKit Mol objects.
    """
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                atom.ClearProp('molAtomMapNumber')
