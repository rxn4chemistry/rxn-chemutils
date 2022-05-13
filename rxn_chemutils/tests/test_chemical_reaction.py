# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2021
# ALL RIGHTS RESERVED

import pytest

from rxn_chemutils.chemical_reaction import ChemicalReaction, ChemicalReactionPart


@pytest.fixture
def chemical_reaction():
    return ChemicalReaction("[14C]Cl.[Na]O>O>[Na]Cl.[14C]O")


@pytest.fixture
def duplicate_chemical_reaction():
    return ChemicalReaction(
        "[14C]Cl.[14C]Cl.[Na]O>O.O>[Na]Cl.[Na]Cl.[14C]O", remove_duplicates=True
    )


@pytest.fixture
def dirty_chemical_reaction():
    return ChemicalReaction("[14C]Cl.[Na]O>O>[Na]Cl.[14C]O.O.[14C]Cl")


def test_remove_duplicates(chemical_reaction, duplicate_chemical_reaction):
    assert chemical_reaction == duplicate_chemical_reaction


def test_len(chemical_reaction):
    assert len(chemical_reaction) == 5


def test_str(chemical_reaction):
    assert str(chemical_reaction) == "[14C]Cl.O[Na]>O>[Na]Cl.[14C]O"


def test_eq(chemical_reaction):
    chemical_reaction_reversed = ChemicalReaction("[Na]Cl.[14C]O>O>[14C]Cl.[Na]O")
    assert chemical_reaction == chemical_reaction
    assert chemical_reaction != chemical_reaction_reversed


def test_reactant_count(chemical_reaction):
    assert len(chemical_reaction.reactants) == 2


def test_agent_count(chemical_reaction):
    assert len(chemical_reaction.agents) == 1


def test_product_count(chemical_reaction):
    assert len(chemical_reaction.products) == 2


def test_get_reactants(chemical_reaction):
    assert chemical_reaction.get_reactants_as_smiles() == ["[14C]Cl", "O[Na]"]


def test_get_agents(chemical_reaction):
    assert chemical_reaction.get_agents_as_smiles() == ["O"]


def test_get_products(chemical_reaction):
    assert chemical_reaction.get_products_as_smiles() == ["[Na]Cl", "[14C]O"]


def test_find(chemical_reaction):
    assert chemical_reaction.find("O") == ([1], [0], [1])


def test_find_in(chemical_reaction):
    assert chemical_reaction.find_in("O", ChemicalReactionPart.reactants) == [1]
    assert chemical_reaction.find_in("O", ChemicalReactionPart.agents) == [0]
    assert chemical_reaction.find_in("O", ChemicalReactionPart.products) == [1]


def test_remove(chemical_reaction):
    chemical_reaction.remove(([1], [0], [1]))
    assert str(chemical_reaction) == "[14C]Cl>>[Na]Cl"


def test_filter(chemical_reaction):
    chemical_reaction.filter(([1], [0], [1]))
    assert str(chemical_reaction) == "O[Na]>O>[14C]O"


def test_sort(chemical_reaction):
    chemical_reaction.sort()
    assert str(chemical_reaction) == "O[Na].[14C]Cl>O>[14C]O.[Na]Cl"


def test_sort_only_reactants(chemical_reaction):
    chemical_reaction.sort(sort_products=False, sort_agents=False)
    assert str(chemical_reaction) == "O[Na].[14C]Cl>O>[Na]Cl.[14C]O"


def test_remove_precursors_from_products(dirty_chemical_reaction):
    dirty_chemical_reaction.remove_precursors_from_products()
    assert str(dirty_chemical_reaction) == "[14C]Cl.O[Na]>O>[Na]Cl.[14C]O"


def test_has_none(chemical_reaction):
    assert not chemical_reaction.has_none()
    chemical_reaction.products.append(None)
    assert chemical_reaction.has_none()


def test_remove_none(chemical_reaction):
    chemical_reaction.products.append(None)
    chemical_reaction.remove_none()
    assert len(chemical_reaction.products) == 2
