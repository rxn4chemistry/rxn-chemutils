from rxn_chemutils.conversion import mol_to_smiles, smiles_to_mol
from rxn_chemutils.rdkit_utils import combine_mols


def test_combine_mols():
    # No molecule -> empty Mol
    assert mol_to_smiles(combine_mols([])) == ""

    # Test a few single mols
    for smiles in ["CC", "CCO", "CC(C)CCO"]:
        mol = smiles_to_mol(smiles)
        assert mol_to_smiles(combine_mols([mol])) == smiles

    # Test a few combinations of molecules
    for multi_smiles in [
        ("CC", "O"),
        ("CC", "CCO", "N"),
        ("CC", "CC", "CC(C)CCO", "c1ccccc1"),
    ]:
        mols = (smiles_to_mol(smiles) for smiles in multi_smiles)
        merged_mol = combine_mols(mols)
        # we verify the merge by checking that the conversion back to SMILES
        # leads to the '.'-concatenated separate SMILES.
        assert mol_to_smiles(merged_mol) == ".".join(multi_smiles)
