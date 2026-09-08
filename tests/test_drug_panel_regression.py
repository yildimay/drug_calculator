"""
Regression tests against the curated real-drug panel (tests/drug_panel.py).

Unlike the single hand-picked examples in test_molecule_parser.py/
test_decision_tree.py/test_admet.py, this runs every branch of the pipeline
against ~17 real, structurally diverse drugs in one sweep, so a change that
breaks one drug class (e.g. heavy-element handling, or the large-non-aromatic
branch) is caught even if the more targeted unit tests happen not to cover it.
"""

import pytest

from core.admet import analyze_drug_likeness
from core.molecule_parser import parse_smiles
from dft.decision_tree import determine_method_and_basis
from tests.drug_panel import DRUG_PANEL


@pytest.mark.parametrize("entry", DRUG_PANEL, ids=[e.name for e in DRUG_PANEL])
def test_panel_entry_parses_successfully(entry):
    molecule_data = parse_smiles(entry.smiles)
    assert molecule_data is not None
    assert molecule_data["valid"] is True, molecule_data.get("error")


@pytest.mark.parametrize("entry", DRUG_PANEL, ids=[e.name for e in DRUG_PANEL])
def test_panel_entry_molecular_properties_match(entry):
    molecule_data = parse_smiles(entry.smiles)

    assert molecule_data["has_transition_metals"] == entry.expected_has_transition_metals
    assert molecule_data["has_aromatic_rings"] == entry.expected_has_aromatic_rings
    assert molecule_data["has_heavy_elements"] == entry.expected_has_heavy_elements
    assert molecule_data["formal_charge"] == entry.expected_formal_charge


@pytest.mark.parametrize("entry", DRUG_PANEL, ids=[e.name for e in DRUG_PANEL])
def test_panel_entry_dft_recommendation_matches(entry):
    molecule_data = parse_smiles(entry.smiles)
    dft_params = determine_method_and_basis(molecule_data)

    assert "error" not in dft_params
    assert dft_params["functional"] == entry.expected_functional
    assert dft_params["basis_set"] == entry.expected_basis_set


@pytest.mark.parametrize("entry", DRUG_PANEL, ids=[e.name for e in DRUG_PANEL])
def test_panel_entry_drug_likeness_matches(entry):
    admet = analyze_drug_likeness(entry.smiles)

    assert admet is not None
    assert admet["lipinski_pass"] == entry.expected_lipinski_pass
    assert admet["veber_pass"] == entry.expected_veber_pass


def test_panel_covers_all_four_functional_branches():
    functionals = {entry.expected_functional for entry in DRUG_PANEL}
    assert functionals == {"TPSSh", "wB97XD", "PBE0", "B3LYP-D3"}


def test_panel_covers_all_four_basis_set_branches():
    basis_sets = {entry.expected_basis_set for entry in DRUG_PANEL}
    assert basis_sets == {"LANL2DZ", "6-31+G(d,p)", "6-311G(d,p)", "6-31G(d)"}


def test_panel_covers_both_lipinski_outcomes():
    outcomes = {entry.expected_lipinski_pass for entry in DRUG_PANEL}
    assert outcomes == {True, False}


def test_panel_covers_both_veber_outcomes():
    outcomes = {entry.expected_veber_pass for entry in DRUG_PANEL}
    assert outcomes == {True, False}
