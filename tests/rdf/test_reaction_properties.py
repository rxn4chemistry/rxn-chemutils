from rxn.chemutils.rdf.reaction_properties import (
    ReactionProperties,
    find_compounds,
    find_compounds_with_category,
)


def test_instantiate_reaction_properties() -> None:
    # Basic example with a reaction containing a mix of values and list properties
    meta = {
        "RXN:TEST_FILE": "1",
        "RXN:TEST_RXN_CONDITIONS": "2",
        "RXN:TEST_RXN_SEQUENCE_NUM": "3",
        "RXN:TEST_TABLEROW_NUM": "4",
        "RXN:SOME_LIST_PROPERTY(1):BROAD": "5",
        "RXN:SOME_LIST_PROPERTY(1):MEDIUM": "6",
        "RXN:SOME_LIST_PROPERTY(1):NARROW": "7",
        "RXN:SOME_LIST_PROPERTY(2):BROAD": "8",
        "RXN:SOME_LIST_PROPERTY(2):MEDIUM": "9",
        "RXN:SOME_LIST_PROPERTY(2):NARROW": "10",
        "RXN:PRODUCT(1):YIELD_%": "11",
        "RXN:PRODUCT(2):YIELD_%": "12",
        "RXN:PRODUCT(4):YIELD_%": "13",
    }
    properties = ReactionProperties(meta)
    assert properties.properties == {
        "TEST_FILE": "1",
        "TEST_RXN_CONDITIONS": "2",
        "TEST_RXN_SEQUENCE_NUM": "3",
        "TEST_TABLEROW_NUM": "4",
        "SOME_LIST_PROPERTY": [
            {
                "BROAD": "5",
                "MEDIUM": "6",
                "NARROW": "7",
            },
            {
                "BROAD": "8",
                "MEDIUM": "9",
                "NARROW": "10",
            },
        ],
        "PRODUCT": [
            {"YIELD_%": "11"},
            {"YIELD_%": "12"},
            {},
            {"YIELD_%": "13"},
        ],
    }


def test_generic_properties() -> None:
    # Example containing some properties from the generic reaction as well
    meta = {
        "RXN:P1": "1",
        "RXN:SOME_P": "2",
        "RXN:SEQ_P": "3",
        "RXN:ROW": "4",
        "RXN:GENERIC:DF": "5",
        "RXN:GENERIC:DUMMY": "6",
    }
    properties = ReactionProperties(meta)
    assert properties.properties == {
        "P1": "1",
        "SOME_P": "2",
        "SEQ_P": "3",
        "ROW": "4",
        "GENERIC": {
            "DF": "5",
            "DUMMY": "6",
        },
    }


def test_prunes_mol_and_symbol() -> None:
    # Example removing "MOL(1):" and "SYMBOL(1):"
    meta = {
        "RXN:SOLVENT(1):MOL(1):MOLSTRUCTURE": "1",
        "RXN:SOLVENT(1):MOL(1):SYMBOL(1):SYMBOL": "2",
        "RXN:SOLVENT(2):MOL(1):MOLSTRUCTURE": "3",
        "RXN:SOLVENT(2):MOL(1):SYMBOL(1):SYMBOL": "4",
        "RXN:GENERIC:SOLVENT(1):MOL(1):MOLSTRUCTURE": "5",
        "RXN:GENERIC:SOLVENT(1):MOL(1):SYMBOL(1):SYMBOL": "6",
    }
    properties = ReactionProperties(meta)
    assert properties.properties == {
        "SOLVENT": [
            {
                "MOLSTRUCTURE": "1",
                "SYMBOL": "2",
            },
            {
                "MOLSTRUCTURE": "3",
                "SYMBOL": "4",
            },
        ],
        "GENERIC": {
            "SOLVENT": [
                {
                    "MOLSTRUCTURE": "5",
                    "SYMBOL": "6",
                },
            ],
        },
    }


def test_find_compounds_simple() -> None:
    property_dict = {
        "MOLSTRUCTURE": "3",
        "SYMBOL": "4",
    }

    assert list(find_compounds(property_dict)) == [
        {
            "MOLSTRUCTURE": "3",
            "SYMBOL": "4",
        },
    ]


def test_find_compounds_nested() -> None:
    property_dict = {
        "SOLVENT": [
            {
                "MOLSTRUCTURE": "1",
                "SYMBOL": "2",
            },
            {
                "MOLSTRUCTURE": "3",
                "SYMBOL": "4",
            },
        ],
        "GENERIC": {
            "CATALYST": [
                {
                    "MOLSTRUCTURE": "5",
                    "SYMBOL": "6",
                },
            ],
        },
    }

    assert list(find_compounds(property_dict)) == [
        {
            "MOLSTRUCTURE": "1",
            "SYMBOL": "2",
        },
        {
            "MOLSTRUCTURE": "3",
            "SYMBOL": "4",
        },
        {
            "MOLSTRUCTURE": "5",
            "SYMBOL": "6",
        },
    ]


def test_find_compounds_with_category() -> None:
    property_dict = {
        "SOLVENT": [
            {
                "MOLSTRUCTURE": "1",
                "SYMBOL": "2",
            },
            {
                "MOLSTRUCTURE": "3",
                "SYMBOL": "4",
            },
        ],
        "GENERIC": {
            "CATALYST": [
                {
                    "MOLSTRUCTURE": "5",
                    "SYMBOL": "6",
                },
            ],
        },
        "NO_CATEGORY": {  # This one is not in a list -> empty category.
            "MOLSTRUCTURE": "8",
            "SYMBOL": "9",
        },
    }

    assert list(find_compounds_with_category(property_dict)) == [
        (
            {
                "MOLSTRUCTURE": "1",
                "SYMBOL": "2",
            },
            "SOLVENT",
        ),
        (
            {
                "MOLSTRUCTURE": "3",
                "SYMBOL": "4",
            },
            "SOLVENT",
        ),
        (
            {
                "MOLSTRUCTURE": "5",
                "SYMBOL": "6",
            },
            "CATALYST",
        ),
        (
            {
                "MOLSTRUCTURE": "8",
                "SYMBOL": "9",
            },
            "",
        ),
    ]
