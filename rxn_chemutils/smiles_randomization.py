import random

from rdkit import Chem

from rxn_chemutils.conversion import smiles_to_mol, mol_to_smiles


def randomize_smiles_rotated(smiles: str, with_order_reversal: bool = True) -> str:
    """
    Randomize a SMILES string by doing a cyclic rotation of the atomic indices.

    Adapted from https://github.com/GLambard/SMILES-X/blob/758478663030580a363a9ee61c11f6d6448e18a1/SMILESX/augm.py#L19.

    Raises:
        InvalidSmiles: for invalid molecules.

    Args:
        smiles: SMILES string to randomize.
        with_order_reversal: whether to reverse the atom order with 50% chance.

    Returns:
        Randomized SMILES string.
    """

    mol = smiles_to_mol(smiles, sanitize=False)

    n_atoms = mol.GetNumAtoms()

    # Generate random values
    rotation_index = random.randint(0, n_atoms - 1)
    reverse_order = with_order_reversal and random.choice([True, False])

    # Generate new atom indices order
    atoms = list(range(n_atoms))
    new_atoms_order = (atoms[rotation_index % len(atoms):] + atoms[:rotation_index % len(atoms)])
    if reverse_order:
        new_atoms_order.reverse()

    mol = Chem.RenumberAtoms(mol, new_atoms_order)
    return mol_to_smiles(mol, canonical=False)


def randomize_smiles_restricted(smiles: str) -> str:
    """
    Randomize a SMILES string in a restricted fashion.

    Raises:
        InvalidSmiles: for invalid molecules.

    Args:
        smiles: SMILES string to randomize.

    Returns:
        Randomized SMILES string.
    """
    mol = smiles_to_mol(smiles, sanitize=False)
    new_atom_order = list(range(mol.GetNumAtoms()))
    random.shuffle(new_atom_order)
    mol = Chem.RenumberAtoms(mol, newOrder=new_atom_order)
    return mol_to_smiles(mol, canonical=False)


def randomize_smiles_unrestricted(smiles: str) -> str:
    """
    Randomize a SMILES string in an unrestricted fashion.

    Raises:
        InvalidSmiles: for invalid molecules.

    Args:
        smiles: SMILES string to randomize.

    Returns:
        Randomized SMILES string.
    """
    mol = smiles_to_mol(smiles, sanitize=False)
    return Chem.MolToSmiles(mol, canonical=False, doRandom=True)
