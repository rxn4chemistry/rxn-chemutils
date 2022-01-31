# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

import pytest
from rdkit import Chem

from rxn_chemutils.conversion import (
    smiles_to_mol, inchify_smiles, canonicalize_smiles, canonicalize_smiles_with_fragment_bonds,
    split_smiles_and_fragment_info, cleanup_smiles, sanitize_mol, mol_to_smiles,
    maybe_canonicalize, smiles_to_inchi, remove_hydrogens, mol_to_mdl, mdl_to_mol
)
from rxn_chemutils.exceptions import InvalidSmiles, SanitizationError, InvalidMdl

WATER_MDL = '''
     RDKit          2D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
'''

CFC_MDL = '''
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
'''


def test_smiles_to_mol():
    mol = smiles_to_mol('COC')
    assert mol.GetNumAtoms() == 3
    assert mol.GetAtomWithIdx(0).GetSymbol() == 'C'
    assert mol.GetAtomWithIdx(1).GetSymbol() == 'O'
    assert mol.GetAtomWithIdx(2).GetSymbol() == 'C'

    # must raise for invalid SMILES (RDKit's default: return None)
    with pytest.raises(InvalidSmiles):
        smiles_to_mol('')

    # Invalid SMILES symbol
    with pytest.raises(InvalidSmiles):
        smiles_to_mol('L')

    # Invalid valence
    with pytest.raises(InvalidSmiles):
        smiles_to_mol('CFC')

    # Radicals are kept when converting back to a SMILES
    assert mol_to_smiles(smiles_to_mol('CC[C]CC')) == 'CC[C]CC'


def test_smiles_to_mol_without_sanitization():
    # With sanitization, this would fail
    mol = smiles_to_mol('CFC', sanitize=False)
    assert mol.GetNumAtoms() == 3
    assert mol.GetAtomWithIdx(0).GetSymbol() == 'C'
    assert mol.GetAtomWithIdx(1).GetSymbol() == 'F'
    assert mol.GetAtomWithIdx(2).GetSymbol() == 'C'

    # empty molecule
    with pytest.raises(InvalidSmiles):
        smiles_to_mol('', sanitize=False)

    # invalid SMILES symbol
    with pytest.raises(InvalidSmiles):
        smiles_to_mol('L', sanitize=False)

    # Doing strictly no sanitization with the RDKit functions messes up the
    # radicals - smiles_to_mol takes care of this
    assert mol_to_smiles(smiles_to_mol('CC[C]CC', sanitize=False)) == 'CC[C]CC'


def test_mdl_to_mol():
    # NB: we assert that the molecules are correct with their SMILES notation
    assert mol_to_smiles(mdl_to_mol(WATER_MDL)) == 'O'

    # Empty MDL
    with pytest.raises(InvalidMdl):
        _ = mdl_to_mol('')

    # Invalid valence - raises by default
    with pytest.raises(InvalidMdl):
        _ = mdl_to_mol(CFC_MDL)

    # Invalid valence - does not raise if no sanitization
    assert mol_to_smiles(mdl_to_mol(CFC_MDL, sanitize=False)) == 'CFC'

    # Back-conversion leads to original MDL
    assert mol_to_mdl(mdl_to_mol(WATER_MDL)) == WATER_MDL
    assert mol_to_mdl(mdl_to_mol(CFC_MDL, sanitize=False)) == CFC_MDL


def test_sanitize_mol():

    def example_mol():
        return smiles_to_mol('C1=CC=CC=C1N(=O)=O', sanitize=False)

    # Example 1: do only aromatization
    mol = example_mol()
    sanitize_mol(mol, include_sanitizations=[Chem.SANITIZE_SETAROMATICITY])
    assert mol_to_smiles(mol, canonical=False) == 'c1ccccc1N(=O)=O'

    # Example 2: do only cleanup
    mol = example_mol()
    sanitize_mol(mol, include_sanitizations=[Chem.SANITIZE_CLEANUP])
    assert mol_to_smiles(mol, canonical=False) == 'C1=CC=CC=C1[N+]([O-])=O'

    # Example 3: do both cleanup and aromatization
    mol = example_mol()
    sanitize_mol(mol, include_sanitizations=[Chem.SANITIZE_CLEANUP, Chem.SANITIZE_SETAROMATICITY])
    assert mol_to_smiles(mol, canonical=False) == 'c1ccccc1[N+]([O-])=O'

    # Example 6: do all except aromatization
    mol = example_mol()
    sanitize_mol(mol, exclude_sanitizations=[Chem.SANITIZE_SETAROMATICITY])
    assert mol_to_smiles(mol, canonical=False) == 'C1=CC=CC=C1[N+]([O-])=O'

    # Example 5: do all except cleanup (and properties -> invalid valence otherwise)
    mol = example_mol()
    sanitize_mol(mol, exclude_sanitizations=[Chem.SANITIZE_CLEANUP, Chem.SANITIZE_PROPERTIES])
    assert mol_to_smiles(mol, canonical=False) == 'c1ccccc1N(=O)=O'

    # Raises exception for invalid sanitization
    mol = smiles_to_mol('CFC', sanitize=False)
    with pytest.raises(SanitizationError):
        sanitize_mol(mol, include_sanitizations=[Chem.SANITIZE_PROPERTIES])

    # Not providing any options does the full sanitization
    mol_1 = example_mol()
    mol_2 = example_mol()
    sanitize_mol(mol_1)
    sanitize_mol(mol_2, include_sanitizations=[Chem.SANITIZE_ALL])
    assert mol_to_smiles(mol_1) == mol_to_smiles(mol_2)

    # Cannot specify all included and excluded sanitizations at the same time
    with pytest.raises(ValueError):
        sanitize_mol(
            example_mol(),
            include_sanitizations=[Chem.SANITIZE_ALL],
            exclude_sanitizations=[Chem.SANITIZE_NONE]
        )


def test_remove_hydrogens():
    smiles = '[H]C([H])([H])Oc1ccc2cccc(C([H])([H])C#N)c2c1'
    smiles_without_h = 'COc1ccc2cccc(CC#N)c2c1'

    # Removing the unnecessary hydrogens
    mol = smiles_to_mol(smiles, sanitize=False)
    assert mol_to_smiles(mol, canonical=False) == smiles
    mol = remove_hydrogens(mol)
    assert mol_to_smiles(mol, canonical=False) == smiles_without_h

    # check that it does not sanitize the molecules
    smiles = '[H]C([H])([H])CN(=O)=O'
    smiles_without_h = 'CCN(=O)=O'
    mol = smiles_to_mol(smiles, sanitize=False)
    mol = remove_hydrogens(mol)
    assert mol_to_smiles(mol, canonical=False) == smiles_without_h


def test_inchify_smiles():
    assert inchify_smiles('Clc1nc2sccc2c(=NCc2ccccc2)[nH]1'
                          ) == inchify_smiles('Clc1nc(NCc2ccccc2)c2ccsc2n1')

    # An empty SMILES string must raise an error
    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('')

    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('L')


def test_canonicalize_smiles_with_fragment_bonds():
    canonical_smiles = canonicalize_smiles_with_fragment_bonds('[Mg+2]~O=S(=O)([O-])[O-]')
    assert canonical_smiles == canonicalize_smiles_with_fragment_bonds('O=S(=O)([O-])[O-]~[Mg+2]')


def test_canonicalize():
    assert canonicalize_smiles('CCO') == canonicalize_smiles('OCC')

    # An empty SMILES string must raise an error
    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('')

    # invalid SMILES symbol
    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('L')

    # invalid valence
    with pytest.raises(InvalidSmiles):
        canonicalize_smiles('CFC')

    # Canonicalization is possible for invalid valence if one does not sanitize
    assert canonicalize_smiles('CFC', check_valence=False) == 'CFC'
    assert canonicalize_smiles('F(C)C', check_valence=False) == 'CFC'
    assert canonicalize_smiles('F(C)(C)', check_valence=False) == 'CFC'

    # For invalid valences, still canonicalizes aromaticity, etc.
    assert (
        canonicalize_smiles('C1=CC=CC=C1CFC', check_valence=False) ==
        canonicalize_smiles('c1ccccc1CFC', check_valence=False)
    )
    # For special aromatic cases, still canonicalizes aromaticity, etc.
    assert (
        canonicalize_smiles('[se]1cccc1', check_valence=False) ==
        canonicalize_smiles('[se]1C=CC=C1', check_valence=False)
    )


def test_canonicalize_rotation():
    # In an earlier version, the following SMILES kept switching between two
    # SMILES variants when canonicalizing (back-and-forth between aromatic and
    # non-aromatic.
    smiles = 'c1cc[se]c1'
    canonical_1 = canonicalize_smiles(smiles)
    canonical_2 = canonicalize_smiles(canonical_1)
    assert canonical_2 == canonical_1


def test_canonicalize_removes_unnecessary_hydrogens():
    # in an earlier version, canonicalize_smiles(x) did not behave the same as
    # Chem.MolToSmiles(Chem.MolFromSmiles(x)) when there were implicit hydrogens.
    smiles = '[H]C([H])([H])Oc1ccc2cccc(C([H])([H])C#N)c2c1'

    assert canonicalize_smiles(smiles) == 'COc1ccc2cccc(CC#N)c2c1'


def test_maybe_canonicalize():
    # basic example
    assert maybe_canonicalize('C(C)C') == 'CCC'

    # invalid -> does not change the SMILES
    assert maybe_canonicalize('C(CC') == 'C(CC'
    assert maybe_canonicalize('F(C)C') == 'F(C)C'

    # invalid valence, but flag given -> does the canonicalization
    assert maybe_canonicalize('F(C)C', check_valence=False) == 'CFC'


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
    assert cleanup_smiles('[se]1cccc1') == '[se]1cccc1'
    assert cleanup_smiles('[se]1C=CC=C1') == '[se]1C=CC=C1'

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


def test_smiles_to_inchi():
    # different SMILES representations for phenol (aromatic, non-aromatic)
    # should be converted to the same InChI
    phenol_inchi = 'InChI=1S/C6H6O/c7-6-4-2-1-3-5-6/h1-5,7H'
    for smiles in ['c1ccccc1O', 'C1=CC=CC=C1O']:
        assert smiles_to_inchi(smiles) == phenol_inchi

    # different SMILES representations for 4(1H)-pyrimidinone (tautomers)
    # should be converted to the same InChI
    expected_inchi = 'InChI=1S/C4H4N2O/c7-4-1-2-5-3-6-4/h1-3H,(H,5,6,7)'
    for smiles in ['C1=CN=CNC1=O', 'C1=CNC=NC1=O', 'C1=CN=CN=C1O']:
        assert smiles_to_inchi(smiles) == expected_inchi

    with pytest.raises(InvalidSmiles):
        smiles_to_inchi('')

    with pytest.raises(InvalidSmiles):
        smiles_to_inchi('invalid')

    # Does not raise for invalid valence
    assert smiles_to_inchi('CFC') == 'InChI=1S/C2H6F/c1-3-2/h1-2H3'

    # SMILES that are aromatized or not lead to the same InChi
    assert smiles_to_inchi('C1=CC=CC=C1') == smiles_to_inchi('c1ccccc1')

    # Stereochemistry information is kept: the InChi should be different if the
    # SMILES have different double-bond or R/S stereochemistries
    assert len({smiles_to_inchi(s) for s in ['CC=CC', 'C/C=C/C', r'C/C=C\C']}) == 3
    assert len({smiles_to_inchi(s) for s in ['C(O)(N)C', '[C@H](O)(N)C', '[C@@H](O)(N)C']}) == 3

    # extended tautomer check - where the default options do not detect
    # tautomerism in a 1,3-diketone
    assert smiles_to_inchi('C1C(=O)CC(=O)CC1') != smiles_to_inchi('C1C(=O)C=C(O)CC1')
    assert smiles_to_inchi('C1C(=O)CC(=O)CC1', True) == smiles_to_inchi('C1C(=O)C=C(O)CC1', True)
