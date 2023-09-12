from rxn.chemutils.rdf.property_parser import (
    PropertyParser,
    parse_properties,
    serialize_properties,
)


def test_parse_property_simple_property() -> None:
    parser = PropertyParser()
    parser.parse_property("PROP", "1")
    assert parser.result == {"PROP": "1"}


def test_parse_property_nested_values() -> None:
    parser = PropertyParser()
    parser.parse_property("A:B", "1")
    parser.parse_property("A:C", "2")
    assert parser.result == {
        "A": {"B": "1", "C": "2"},
    }


def test_parse_property_mixed_nested_values() -> None:
    parser = PropertyParser()
    parser.parse_property("A:B", "1")
    parser.parse_property("A:C", "2")
    parser.parse_property("A:D:E", "3")
    parser.parse_property("A:D:F", "4")
    assert parser.result == {
        "A": {
            "B": "1",
            "C": "2",
            "D": {
                "E": "3",
                "F": "4",
            },
        }
    }


def test_parse_property_list() -> None:
    parser = PropertyParser()
    parser.parse_property("A(1):B", "1")
    parser.parse_property("A(1):C", "2")
    parser.parse_property("A(2):B", "3")
    parser.parse_property("A(2):C", "4")
    assert parser.result == {
        "A": [
            {
                "B": "1",
                "C": "2",
            },
            {
                "B": "3",
                "C": "4",
            },
        ]
    }


def test_parse_property_list_incorrect_order() -> None:
    # Adding elements to list not in the right order, and with not all of them specified
    parser = PropertyParser()
    parser.parse_property("A(4):B", "3")
    parser.parse_property("A(4):C", "4")
    parser.parse_property("A(1):B", "1")
    parser.parse_property("A(1):C", "2")
    assert parser.result == {
        "A": [
            {
                "B": "1",
                "C": "2",
            },
            {},
            {},
            {
                "B": "3",
                "C": "4",
            },
        ]
    }


def test_parse_dict() -> None:
    parser = PropertyParser()
    parser.parse_dict(
        {
            "A(1):B": "1",
            "A(1):C": "2",
            "A(2):B": "3",
            "A(2):C": "4",
            "F:G:H": "5",
        }
    )
    assert parser.result == {
        "A": [
            {
                "B": "1",
                "C": "2",
            },
            {
                "B": "3",
                "C": "4",
            },
        ],
        "F": {"G": {"H": "5"}},
    }


def test_roundtrip() -> None:
    # 1) basic property dictionary
    properties = {
        "A(1):B": "1",
        "A(1):C": "2",
        "A(2):B": "3",
        "A(2):C": "4",
        "F:G:H": "5",
    }
    assert serialize_properties(parse_properties(properties)) == properties

    # 2) with incomplete lists (indexes jumping)
    properties = {
        "A(1):B": "1",
        "A(1):C": "2",
        "A(3):B": "3",
        "A(6):C": "4",
        "F:G:H": "5",
    }
    assert serialize_properties(parse_properties(properties)) == properties
