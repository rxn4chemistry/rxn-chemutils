import random

from rdkit import Chem

from .conversion import mol_to_smiles, smiles_to_mol

# Highest value to give as a random seed for RDKit.
# Any value higher than that will cause problems.
_MAX_RDKIT_RANDOM_SEED = 2147483647


def randomize_smiles_rotated(smiles: str, with_order_reversal: bool = True) -> str:
    """
    Randomize a SMILES string by doing a cyclic rotation of the atomic indices.

    Adapted from https://github.com/GLambard/SMILES-X/blob/758478663030580a363a9ee61c11f6d6448e18a1/SMILESX/augm.py#L19.

    The outputs of this function can be reproduced by setting the seed with random.seed().

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
    new_atoms_order = (
        atoms[rotation_index % len(atoms) :] + atoms[: rotation_index % len(atoms)]
    )
    if reverse_order:
        new_atoms_order.reverse()

    mol = Chem.RenumberAtoms(mol, new_atoms_order)
    return mol_to_smiles(mol, canonical=False)


def randomize_smiles_restricted(smiles: str) -> str:
    """
    Randomize a SMILES string in a restricted fashion.

    The outputs of this function can be reproduced by setting the seed with random.seed().

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

    The outputs of this function can be reproduced by setting the seed with random.seed().

    Raises:
        InvalidSmiles: for invalid molecules.

    Args:
        smiles: SMILES string to randomize.

    Returns:
        Randomized SMILES string.
    """
    mol = smiles_to_mol(smiles, sanitize=False)

    # We sample the seed to give to RDKit. This makes the call reproducible
    # if one sets random.seed() outside this function.
    seed = random.randint(1, _MAX_RDKIT_RANDOM_SEED)

    # Note: to allow for reproducibility, we do not rely on
    #       Chem.MolToSmiles(mol, canonical=False, doRandom=True)
    # See https://www.rdkit.org/docs/Cookbook.html#enumerate-smiles
    return Chem.MolToRandomSmilesVect(mol, 1, seed)[0]
