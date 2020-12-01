from typing import List

from rdkit.Chem import Mol

from .conversion import smiles_to_mol


class ChemicalReaction:
    """
    RDKit's ChemicalReaction sometimes has errors when parsing reaction SMILES.
    See https://github.com/rdkit/rdkit/issues/2196 for instance.
    This classes follows its interface, without those errors.
    """

    def __init__(self, reaction_smiles: str, sanitize: bool = True):
        self.reaction_smiles = reaction_smiles
        smiles_groups = reaction_smiles.split('>')
        assert len(smiles_groups) == 3

        mol_groups = [
            self.convert_to_mols(smiles_group, sanitize=sanitize) for smiles_group in smiles_groups
        ]
        self.reactants, self.agents, self.products = mol_groups

    def convert_to_mols(self, smiles_group: str, sanitize: bool = True) -> List[Mol]:
        """
        Args:
            smiles_group: List of SMILES strings, separated by dots
            sanitize: whether to apply RDKit sanitization during conversion
        """
        # Handle empty strings:
        if not smiles_group:
            return []

        smiles_list = smiles_group.split('.')
        return [smiles_to_mol(smiles, sanitize=sanitize) for smiles in smiles_list]

    def GetReactants(self) -> List[Mol]:
        return self.reactants

    def GetAgents(self) -> List[Mol]:
        return self.agents

    def GetProducts(self) -> List[Mol]:
        return self.products
