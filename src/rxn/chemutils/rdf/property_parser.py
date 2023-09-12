import re
from typing import Any, Dict, Generator, List, Tuple

from rxn.utilities.regex import capturing

_LIST_REGEX_STRING = capturing(".+") + r"\(" + capturing(r"\d+") + r"\)"
_LIST_REGEX = re.compile(_LIST_REGEX_STRING)


class PropertyParser:
    """
    To parse properties given in a RDF / MDL into a nested dictionary.
    """

    def __init__(self) -> None:
        self.result: Dict[str, Any] = {}

    def parse_dict(self, property_dict: Dict[str, str]) -> None:
        for key, value in property_dict.items():
            self.parse_property(key, value)

    def parse_property(self, key: str, value: str) -> None:
        self._parse_property(self.result, key, value)

    def _is_list_property(self, subkey: str) -> bool:
        return _LIST_REGEX.match(subkey) is not None

    def _parse_property(self, container: Dict[str, Any], key: str, value: str) -> None:
        if not key:
            raise ValueError("A key must be provided.")

        splits = key.split(":")
        if len(splits) == 1:
            container[key] = value
            return

        if self._is_list_property(splits[0]):
            self._parse_list_property(container, splits, value)
        else:
            self._parse_dict_property(container, splits, value)

    def _parse_list_property(
        self, container: Dict[str, Any], key_splits: List[str], value: str
    ) -> None:
        list_match = _LIST_REGEX.match(key_splits[0])
        if list_match is None:
            raise RuntimeError(
                "Not a list property - by construction, this should not happen."
            )

        subkey = list_match.group(1)
        list_index = int(list_match.group(2))

        # Initialize list if needed
        if subkey not in container:
            container[subkey] = []

        subkey_list = container[subkey]

        # Add new dicts if necessary
        for _ in range(len(subkey_list), list_index):
            subkey_list.append({})

        self._parse_property(
            subkey_list[list_index - 1], ":".join(key_splits[1:]), value
        )

    def _parse_dict_property(
        self, container: Dict[str, Any], key_splits: List[str], value: str
    ) -> None:
        subkey = key_splits[0]
        if subkey not in container:
            container[subkey] = {}
        self._parse_property(container[subkey], ":".join(key_splits[1:]), value)


class PropertySerializer:
    """Do the reverse operation compared to PropertyParser."""

    def convert_dict(self, container: Dict[str, Any]) -> Dict[str, str]:
        return {
            key: value
            for key, value in self._convert_dict(prefix="", container=container)
        }

    def _convert_dict(
        self, prefix: str, container: Dict[str, Any]
    ) -> Generator[Tuple[str, str], None, None]:
        for key, value in container.items():
            yield from self._convert(prefix=prefix, key=key, current=value)

    def _convert(
        self, prefix: str, key: str, current: Any
    ) -> Generator[Tuple[str, str], None, None]:
        if isinstance(current, str):
            yield f"{prefix}{key}", current
        elif isinstance(current, list):
            for index, v in enumerate(current, 1):
                yield from self._convert_dict(
                    prefix=f"{prefix}{key}({index}):", container=v
                )
        elif isinstance(current, dict):
            yield from self._convert_dict(prefix=f"{prefix}{key}:", container=current)
        else:
            raise RuntimeError(f"Not supported for property serialization: {current}")


def parse_properties(properties: Dict[str, str]) -> Dict[str, Any]:
    """Parse the properties given in the RDF into a nested dictionary."""
    pp = PropertyParser()
    pp.parse_dict(properties)
    return pp.result


def serialize_properties(properties: Dict[str, Any]) -> Dict[str, str]:
    """Do the reverse operation compared to parse_properties."""
    return PropertySerializer().convert_dict(properties)
