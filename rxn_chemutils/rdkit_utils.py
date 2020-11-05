from typing import Sequence

from rdkit.Chem.rdchem import Mol


def remove_atom_mapping(mols: Sequence[Mol]) -> None:
    for mol in mols:
        for atom in mol.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                atom.ClearProp('molAtomMapNumber')
