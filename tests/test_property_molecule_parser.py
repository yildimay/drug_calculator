"""
Property-based tests for core/molecule_parser.py using Hypothesis.

Two kinds of properties are checked here:
1. Robustness: arbitrary text (not necessarily valid SMILES) must never
   crash parse_smiles() - it should always come back as a well-formed
   "invalid" result instead of raising.
2. Correctness invariants over the real, verified drug panel (sampled_from,
   not randomly generated - generating valid random SMILES from scratch is
   its own hard problem and out of scope here).
"""

from functools import lru_cache

from hypothesis import HealthCheck, given, settings, strategies as st

from core.molecule_parser import parse_smiles
from tests.drug_panel import DRUG_PANEL

_PANEL_SMILES = [entry.smiles for entry in DRUG_PANEL]

# 3D embedding/MMFF optimization for the panel's larger molecules (paclitaxel,
# the cyclosporine-like peptide) is not free - cache so each distinct SMILES
# is only embedded once across all the property checks below, rather than
# once per (test function x Hypothesis example).
_cached_parse_smiles = lru_cache(maxsize=None)(parse_smiles)

# sampled_from over a ~17-item fixed pool doesn't benefit from Hypothesis's
# usual hundreds of examples (there's no infinite input space to explore) -
# a handful is enough to exercise the whole panel across shrinking/reruns.
_PANEL_SETTINGS = settings(max_examples=20, deadline=None)


@given(st.text(max_size=50))
@settings(suppress_health_check=[HealthCheck.too_slow], deadline=None)
def test_arbitrary_text_never_crashes_parse_smiles(text):
    result = parse_smiles(text)
    assert isinstance(result, dict)
    assert "valid" in result
    if not result["valid"]:
        assert "error" in result


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_panel_smiles_always_parse_as_valid(smiles):
    result = _cached_parse_smiles(smiles)
    assert result["valid"] is True


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_valid_result_always_has_required_fields(smiles):
    result = _cached_parse_smiles(smiles)
    required_fields = {
        "smiles", "valid", "formal_charge", "spin_multiplicity", "num_heavy_atoms",
        "has_aromatic_rings", "has_transition_metals", "has_heavy_elements",
        "molecular_weight", "num_atoms_total", "xyz_coordinates", "atom_list",
    }
    assert required_fields.issubset(result.keys())


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_spin_multiplicity_always_positive(smiles):
    result = _cached_parse_smiles(smiles)
    assert result["spin_multiplicity"] >= 1


@given(st.sampled_from(_PANEL_SMILES))
@_PANEL_SETTINGS
def test_num_atoms_total_matches_atom_list_length(smiles):
    result = _cached_parse_smiles(smiles)
    assert result["num_atoms_total"] == len(result["atom_list"])


@given(st.sampled_from(_PANEL_SMILES), st.integers(min_value=-3, max_value=3))
@_PANEL_SETTINGS
def test_charge_override_always_applied_exactly(smiles, charge):
    result = parse_smiles(smiles, charge_override=charge)
    assert result["valid"] is True
    assert result["formal_charge"] == charge


@given(st.sampled_from(_PANEL_SMILES), st.integers(min_value=1, max_value=5))
@_PANEL_SETTINGS
def test_multiplicity_override_always_applied_exactly(smiles, multiplicity):
    result = parse_smiles(smiles, multiplicity_override=multiplicity)
    assert result["valid"] is True
    assert result["spin_multiplicity"] == multiplicity
