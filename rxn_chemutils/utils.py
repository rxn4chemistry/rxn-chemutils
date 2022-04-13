"""Simple utils not involving RDKit."""

# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

import re


def remove_atom_mapping(smiles: str) -> str:
    """
    Remove the atom mapping of a reaction SMILES.

    The resulting SMILES strings will still contain brackets and it may be
    advisable to canonicalize them or clean them up as a postprocessing step.

    Args:
        smiles: SMILES string potentially containing mapping information.

    Returns:
        A SMILES string without atom mapping information.
    """

    # We look for ":" followed by digits before a "]" not coming after an "*"
    return re.sub(r'(?<=[^\*])(:\d+)]', ']', smiles)
