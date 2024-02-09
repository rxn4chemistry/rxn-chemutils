import logging
from typing import Optional

from .conversion import inchi_to_mol, mol_to_inchi, mol_to_smiles, smiles_to_mol
from .exceptions import InvalidInchi, InvalidSmiles

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

RXN_SMILES_SEPARATOR = ">>"


def standardize_smiles(
    smiles: str,
    canonicalize: bool = True,
    sanitize: bool = True,
    find_radicals: bool = True,
    inchify: bool = False,
) -> str:
    """Ensure that a SMILES follows a desired standard.

    It allows canonicalization, sanitization and inchification keeping stereochemistry with isomericSmile=True.
    It can process multiple molecules separated by ".".
    Note that inchify set to True will also canonicalize the molecule.

    Args:
        smiles: SMILES representation of a molecule.
        canonicalize: canonicalize SMILES. Defaults to True.
        sanitize: sanitize SMILES. Defaults to True.
        inchify: inchify the SMILES. Defaults to False.

    Returns:
        a SMILES following the desired standard.
    """
    try:
        molecule = smiles_to_mol(smiles, sanitize=sanitize, find_radicals=find_radicals)
    except InvalidSmiles:
        logger.error(f"SMILES parsing failure: {smiles}.")
        raise

    if inchify:
        try:
            inchi_string = mol_to_inchi(molecule)
        except InvalidInchi:
            logger.error(
                f"Inchification failure for SMILES: {smiles}. Returning its canonical version."
            )
            return mol_to_smiles(molecule, isomericSmiles=True)
        else:
            # canonical set to True because we can't guarantee no canonicalization
            try:
                molecule_from_inchi = inchi_to_mol(inchi_string)
            except InvalidInchi:
                logger.error(
                    f"De-inchification failure for InChi: {inchi_string}. Returning its canonical version."
                )
                return mol_to_smiles(molecule, isomericSmiles=True)
            return mol_to_smiles(molecule_from_inchi, canonical=True)
    if canonicalize:
        return mol_to_smiles(molecule, isomericSmiles=True)
    else:
        return smiles


def standardize_molecules(
    molecules: str,
    canonicalize: bool = True,
    sanitize: bool = True,
    inchify: bool = False,
    fragment_bond: str = "~",
    ordered_precursors: bool = True,
    molecule_token_delimiter: Optional[str] = None,
    is_enzymatic: bool = False,
    enzyme_separator: str = "|",
) -> str:
    """Ensure that a set of molecules represented by a string follows a desired standard.

    Args:
        molecules: molecules SMILES. Molecules can be separated via a ".".
            Fragments are supported with a custom `fragment_bond`.
        canonicalize: canonicalize SMILES. Defaults to True.
        sanitize: sanitize SMILES. Defaults to True.
        inchify: inchify the SMILES. Defaults to False.
        fragment_bond: fragment bond. Defaults to '~'.
        ordered_precursors: order precursors. Defaults to True.
        molecule_token_delimiter: delimiter for big molecule tokens. Defaults to None
        is_enzymatic: the molecules are representing an enzymatic reaction. Defaults to False.
        enzyme_separator: separator for molecules and the enzyme. Defaults to '|'.

    Returns:
        standardized molecules.

    Examples:
        Standardize multiple molecules:
        >>> standardize_molecules('CCO.CC')
        'CC.CCO'
        Standardize multiple molecules including fragment information:
        >>> standardize_molecules('CCO.CC~C')
        'CCO.C~CC'
    """
    enzyme = ""
    if is_enzymatic:
        splitted_molecules = molecules.split(enzyme_separator)
        molecules = splitted_molecules[0]
        if len(splitted_molecules) > 1:
            enzyme = splitted_molecules[1]
            enzyme = "{}{}".format(enzyme_separator, enzyme)
    if molecule_token_delimiter is not None:
        molecules = molecules.replace(molecule_token_delimiter, "")
    if fragment_bond in molecules:
        standardized_molecules_list = [
            # make sure we remove the fragment to have valid SMILES
            standardize_smiles(
                molecule.replace(fragment_bond, "."),
                canonicalize=canonicalize,
                sanitize=sanitize,
                inchify=inchify,
            ).replace(".", fragment_bond)
            for molecule in molecules.split(".")
        ]
        if ordered_precursors:
            standardized_molecules_list = sorted(standardized_molecules_list)
        standardized_molecules = ".".join(standardized_molecules_list)
    else:
        if ordered_precursors:
            # RDKit guarantees ordered precursors
            standardized_molecules = standardize_smiles(
                molecules,
                canonicalize=canonicalize,
                sanitize=sanitize,
                inchify=inchify,
            )
        else:
            standardized_molecules_list = [
                standardize_smiles(
                    molecule,
                    canonicalize=canonicalize,
                    sanitize=sanitize,
                    inchify=inchify,
                )
                for molecule in molecules.split(".")
            ]
            standardized_molecules = ".".join(standardized_molecules_list)
    # add optional enzyme information
    standardized_molecules = "{}{}".format(standardized_molecules, enzyme)
    return standardized_molecules
