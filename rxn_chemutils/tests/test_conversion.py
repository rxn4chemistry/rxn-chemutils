import pytest

from rxn_chemutils.conversion import (
    smiles_to_mol, InvalidSmiles, inchify_smiles, canonicalize_smiles,
    inchify_smiles_with_fragment_bonds, canonicalize_smiles_with_fragment_bonds,
    canonicalize_reaction_smiles, split_smiles_and_fragment_info, cleanup_smiles
)


def test_smiles_to_mol():
    mol = smiles_to_mol('COC')
    assert mol.GetNumAtoms() == 3
    assert mol.GetAtomWithIdx(0).GetSymbol() == 'C'
    assert mol.GetAtomWithIdx(1).GetSymbol() == 'O'
    assert mol.GetAtomWithIdx(2).GetSymbol() == 'C'

    # must raise for invalid SMILES (RDKit's default: return None)
    with pytest.raises(InvalidSmiles):
        smiles_to_mol('')

    with pytest.raises(InvalidSmiles):
        smiles_to_mol('L')


def test_inchify_smiles():
    assert inchify_smiles('Clc1nc2sccc2c(=NCc2ccccc2)[nH]1'
                          ) == inchify_smiles('Clc1nc(NCc2ccccc2)c2ccsc2n1')

    # An empty SMILES string must raise an error
    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('')

    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('L')


def test_inchify_smiles_with_fragment_bonds():
    inchify_smiles = inchify_smiles_with_fragment_bonds('[Mg+2]~O=S(=O)([O-])[O-]')
    assert inchify_smiles == inchify_smiles_with_fragment_bonds('O=S(=O)([O-])[O-]~[Mg+2]')


def test_canonicalize_smiles_with_fragment_bonds():
    canonical_smiles = canonicalize_smiles_with_fragment_bonds('[Mg+2]~O=S(=O)([O-])[O-]')
    assert canonical_smiles == canonicalize_smiles_with_fragment_bonds('O=S(=O)([O-])[O-]~[Mg+2]')


def test_canonicalize():
    assert canonicalize_smiles('CCO') == canonicalize_smiles('OCC')

    # An empty SMILES string must raise an error
    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('')

    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('L')


def test_canonicalize_reaction_smiles():
    reaction_smiles = '[CH:8]1=[CH:7][C:4](=[CH:3][C:2](=[CH:9]1)[Br:1])[CH:5]=[O:6].[CH2:10]([CH2:11]O)[OH:12]' \
                      '>CC1=CC=CC=C1.CC1=CC=C(C=C1)S(=O)(=O)O' \
                      '>[CH2:10]1[CH2:11][O:6][CH:5]([O:12]1)[C:4]2=[CH:3][C:2](=[CH:9][CH:8]=[CH:7]2)[Br:1]'

    canonical = canonicalize_reaction_smiles(reaction_smiles)

    assert canonical == 'O=Cc1cccc(Br)c1.OCCO' \
                        '>Cc1ccccc1.Cc1ccc(S(=O)(=O)O)cc1' \
                        '>Brc1cccc(C2OCCO2)c1'


def test_canonicalize_reaction_smiles_with_fragments():
    # No change to the fragment information should happen
    reaction_smiles = '[CH3:1][C:2]([CH3:3])([CH3:4])[S@@:5](=[O:6])[NH2:7].[CH:14]1=[CH:13][C:11](=[CH:10][C:9](=[CH:15]1)[Br:8])[CH:12]=O' \
                      '>C(Cl)Cl.C(Cl)Cl.C(=O)([O-])[O-].[Cs+].[Cs+]' \
                      '>[CH3:1][C:2]([CH3:3])([CH3:4])[S@@:5](=[O:6])/[N:7]=[CH:12]/[C:11]1=[CH:10][C:9](=[CH:15][CH:14]=[CH:13]1)[Br:8]' \
                      ' |f:4.5.6|'

    canonical = canonicalize_reaction_smiles(reaction_smiles)

    assert canonical == 'CC(C)(C)[S@](N)=O.O=Cc1cccc(Br)c1' \
                        '>ClCCl.ClCCl.O=C([O-])[O-].[Cs+].[Cs+]' \
                        '>CC(C)(C)[S@@](=O)/N=C/c1cccc(Br)c1' \
                        ' |f:4.5.6|'


def test_cleanup_smiles():
    # Basic checks
    assert cleanup_smiles('C') == 'C'
    assert cleanup_smiles('[C]') == '[C]'
    assert cleanup_smiles('[CH4]') == 'C'
    assert cleanup_smiles('C[C]C') == 'C[C]C'
    assert cleanup_smiles('C[CH2]C') == 'CCC'

    # Removes unnecessary dashes and parentheses
    assert cleanup_smiles('C-C') == 'CC'
    assert cleanup_smiles('C(C)') == 'CC'

    # Does not do any canonicalization
    assert cleanup_smiles('C(C)O') == 'C(C)O'

    # Does not do any kekulization / aromatization
    assert cleanup_smiles('c1ccccc1') == 'c1ccccc1'
    assert cleanup_smiles('C1=CC=CC=C1') == 'C1=CC=CC=C1'

    # Does not mess up with fragment molecules
    assert cleanup_smiles('[C].[Pd]') == '[C].[Pd]'
    assert cleanup_smiles('[Pd].[C]') == '[Pd].[C]'

    # Does no valence check
    carbon_with_eight_neighbours = 'C(C)(C)(C)(C)(C)(C)(C)C'
    assert cleanup_smiles(carbon_with_eight_neighbours) == carbon_with_eight_neighbours
    na_with_two_neighbors = '[Na]1OS(=O)(=O)O1'
    assert cleanup_smiles(na_with_two_neighbors) == na_with_two_neighbors

    # More complex example. It is not possible to avoid some reordering here,
    # it should actually be CC(C)c1ccc(cc1)C(=O)CCCCl.
    assert cleanup_smiles(
        '[CH3][CH]([CH3])[c]1[cH][cH][c]([cH][cH]1)[C](=[O])[CH2][CH2][CH2][Cl]'
    ) == 'CC(C)c1ccc(C(=O)CCCCl)cc1'

    # Raises for invalid atom types
    with pytest.raises(InvalidSmiles):
        cleanup_smiles('CC[Invalid]O')

    # Raises for empty string
    with pytest.raises(InvalidSmiles):
        cleanup_smiles('')


def test_split_reaction_smiles():
    reaction_smiles = 'BrC1=C(C(=O)O)C=C(C=C1)COC.C1(O)=CC(O)=CC=C1.[OH-].[Na+]' \
                      '>S(=O)(=O)([O-])[O-].[Cu+2].O' \
                      '>OC1=CC=C2C3=C(C(OC2=C1)=O)C=C(C=C3)COC' \
                      ' |f:2.3,4.5|'

    pure_smiles, fragment_info = split_smiles_and_fragment_info(reaction_smiles)
    assert pure_smiles == 'BrC1=C(C(=O)O)C=C(C=C1)COC.C1(O)=CC(O)=CC=C1.[OH-].[Na+]' \
                          '>S(=O)(=O)([O-])[O-].[Cu+2].O' \
                          '>OC1=CC=C2C3=C(C(OC2=C1)=O)C=C(C=C3)COC'
    assert fragment_info == '|f:2.3,4.5|'


def test_split_reaction_smiles_with_no_fragment_info():
    reaction_smiles = '[CH3:7][C:4]1=[CH:3][C:2](=[C:1]([CH:6]=[CH:5]1)Br)[F:8]' \
                      '>[Li]CCCC.C1CCOC1' \
                      '>[CH3:7][C:4]1=[CH:3][C:2](=[C:1]([CH:6]=[CH:5]1)C(=O)O)[F:8]'

    pure_smiles, fragment_info = split_smiles_and_fragment_info(reaction_smiles)
    assert pure_smiles == reaction_smiles
    assert fragment_info == ''
