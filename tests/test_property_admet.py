"""
Property-based tests for core/admet.py using Hypothesis.
"""

from hypothesis import HealthCheck, given, settings, strategies as st

from core.admet import LIPINSKI_MAX_TOLERATED_VIOLATIONS, analyze_drug_likeness
from tests.drug_panel import DRUG_PANEL

_PANEL_SMILES = [entry.smiles for entry in DRUG_PANEL]
_PANEL_SETTINGS = settings(max_examples=20, deadline=None)


@given(st.text(max_size=50))
@settings(suppress_health_check=[HealthCheck.too_slow], deadline=None)
def test_arbitrary_text_never_crashes_analyze_drug_likeness(text):
    result = analyze_drug_likeness(text)
    assert result is None or isinstance(result, dict)


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_panel_smiles_always_produce_a_result(smiles):
    result = analyze_drug_likeness(smiles)
    assert result is not None


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_descriptor_values_always_non_negative(smiles):
    result = analyze_drug_likeness(smiles)
    assert result["molecular_weight"] >= 0
    assert result["num_h_donors"] >= 0
    assert result["num_h_acceptors"] >= 0
    assert result["tpsa"] >= 0
    assert result["num_rotatable_bonds"] >= 0


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_lipinski_pass_matches_violation_count(smiles):
    # lipinski_pass must always be exactly "violations <= tolerance", not
    # some independently-drifted computation.
    result = analyze_drug_likeness(smiles)
    expected_pass = len(result["lipinski_violations"]) <= LIPINSKI_MAX_TOLERATED_VIOLATIONS
    assert result["lipinski_pass"] == expected_pass


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_veber_pass_matches_violation_count(smiles):
    result = analyze_drug_likeness(smiles)
    expected_pass = len(result["veber_violations"]) == 0
    assert result["veber_pass"] == expected_pass


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_suggestions_present_iff_a_rule_fails(smiles):
    result = analyze_drug_likeness(smiles)
    any_failed = not result["lipinski_pass"] or not result["veber_pass"]
    assert bool(result["suggestions"]) == any_failed


