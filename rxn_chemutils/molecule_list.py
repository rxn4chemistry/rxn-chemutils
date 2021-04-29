from collections import UserList
from typing import TypeVar, Optional, Type, List

T = TypeVar('T', bound='MoleculeList')


class MoleculeList(UserList):
    """
    Container for a set of SMILES string.

    Use this class instead of a simple Python list if the conversion to and from
    SMILES for the whole set is necessary.
    """

    def to_list(self) -> List[str]:
        return self.data

    def to_string(self, fragment_bond: Optional[str] = None) -> str:
        """
        Convert a MoleculeList to a string representation.
        """
        molecules = self.data

        # replace fragment bonds if necessary
        if fragment_bond is not None:
            molecules = [molecule.replace('.', fragment_bond) for molecule in molecules]

        return '.'.join(molecules)

    @classmethod
    def from_string(cls: Type[T], molecules_string: str, fragment_bond: Optional[str] = None) -> T:
        """
        Convert a MoleculeList from a string representation.
        """

        molecules = molecules_string.split('.')
        molecules = [molecule for molecule in molecules if molecule != '']

        # replace fragment bonds if necessary
        if fragment_bond is not None:
            molecules = [molecule.replace(fragment_bond, '.') for molecule in molecules]

        return cls(molecules)


def molecule_list_from_string(molecules_string: str,
                              fragment_bond: Optional[str] = None) -> List[str]:
    """
    Function to convert a string of molecules into a list of molecules (taking
    fragment bonds into account).
    """
    return MoleculeList.from_string(molecules_string, fragment_bond).to_list()


def molecule_list_to_string(molecule_list: List[str], fragment_bond: Optional[str] = None) -> str:
    """
    Function to convert a list of molecules into a string representation
    (taking fragment bonds into account).
    """
    return MoleculeList(molecule_list).to_string(fragment_bond)
