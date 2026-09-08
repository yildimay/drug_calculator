"""
Property-based tests for dft/decision_tree.py using Hypothesis.

Unlike the example-based tests in test_decision_tree.py (a handful of fixed
molecule_data dicts), these generate many random combinations of molecular
property flags and check invariants that must hold for *any* input, not
just the cases we thought to write down by hand. This is the highest-value
place in the codebase for this kind of testing: determine_method_and_basis()
takes a plain dict of primitives (no RDKit/SMILES needed), so Hypothesis can
exercise it directly and exhaustively without any embedding/parsing cost.
"""

from hypothesis import given, strategies as st

from dft.decision_tree import determine_method_and_basis

KNOWN_FUNCTIONALS = {"TPSSh", "wB97XD", "PBE0", "B3LYP-D3"}
KNOWN_BASIS_SETS = {"LANL2DZ", "6-31+G(d,p)", "6-311G(d,p)", "6-31G(d)"}

_ATOM_SYMBOLS = ["C", "H", "O", "N", "F", "Cl", "S", "P"]

molecule_data_strategy = st.fixed_dictionaries({
    "valid": st.just(True),
    "has_transition_metals": st.booleans(),
    "has_aromatic_rings": st.booleans(),
    "has_heavy_elements": st.booleans(),
    "formal_charge": st.integers(min_value=-4, max_value=4),
    "num_atoms_total": st.integers(min_value=1, max_value=400),
    "num_heavy_atoms": st.integers(min_value=1, max_value=250),
    "atom_list": st.lists(
        st.tuples(st.sampled_from(_ATOM_SYMBOLS), st.just(0.0), st.just(0.0), st.just(0.0)),
        min_size=0,
        max_size=40,
    ),
})


@given(molecule_data_strategy)
def test_never_raises_and_returns_no_error_for_valid_input(molecule_data):
    result = determine_method_and_basis(molecule_data)
    assert "error" not in result


@given(molecule_data_strategy)
def test_functional_and_basis_set_always_known(molecule_data):
    result = determine_method_and_basis(molecule_data)
    assert result["functional"] in KNOWN_FUNCTIONALS
    assert result["basis_set"] in KNOWN_BASIS_SETS


@given(molecule_data_strategy)
def test_transition_metals_always_win_functional_selection(molecule_data):
    # Highest-priority rule: transition metals present -> TPSSh, no matter
    # what else is true about the molecule.
    result = determine_method_and_basis(molecule_data)
    if molecule_data["has_transition_metals"]:
        assert result["functional"] == "TPSSh"


@given(molecule_data_strategy)
def test_aromatic_wins_over_size_when_no_transition_metal(molecule_data):
    # Second-priority rule: aromatic (and no transition metal) -> wB97XD,
    # regardless of size - this is the exact behavior that made paclitaxel
    # (62 heavy atoms, aromatic) get wB97XD instead of PBE0.
    result = determine_method_and_basis(molecule_data)
    if not molecule_data["has_transition_metals"] and molecule_data["has_aromatic_rings"]:
        assert result["functional"] == "wB97XD"


@given(molecule_data_strategy)
def test_large_nonaromatic_nonmetal_gets_pbe0(molecule_data):
    result = determine_method_and_basis(molecule_data)
    if (
        not molecule_data["has_transition_metals"]
        and not molecule_data["has_aromatic_rings"]
        and molecule_data["num_heavy_atoms"] > 50
    ):
        assert result["functional"] == "PBE0"


@given(molecule_data_strategy)
def test_small_nonaromatic_nonmetal_gets_default(molecule_data):
    result = determine_method_and_basis(molecule_data)
    if (
        not molecule_data["has_transition_metals"]
        and not molecule_data["has_aromatic_rings"]
        and molecule_data["num_heavy_atoms"] <= 50
    ):
        assert result["functional"] == "B3LYP-D3"


@given(molecule_data_strategy)
def test_heavy_elements_always_win_basis_set_selection(molecule_data):
    result = determine_method_and_basis(molecule_data)
    if molecule_data["has_heavy_elements"]:
        assert result["basis_set"] == "LANL2DZ"
        assert result["ecp_used"] is True


@given(molecule_data_strategy)
def test_hardware_estimate_always_positive(molecule_data):
    result = determine_method_and_basis(molecule_data)
    assert result["nprocshared"] > 0
    assert result["mem_per_proc"] > 0
    assert result["total_mem"] > 0
    assert result["total_mem"] == result["nprocshared"] * result["mem_per_proc"]


@given(molecule_data_strategy)
def test_rationale_is_nonempty_and_mentions_the_picks(molecule_data):
    result = determine_method_and_basis(molecule_data)
    assert result["rationale"]
    assert result["functional"] in result["rationale"]
    assert result["basis_set"] in result["rationale"]


@given(
    molecule_data_strategy,
    st.text(min_size=1, max_size=30),
    st.text(min_size=1, max_size=30),
)
def test_overrides_always_take_effect_verbatim(molecule_data, functional_override, basis_set_override):
    # An override must be used exactly as given, never silently modified or
    # rejected, no matter what the auto-selection would have picked.
    result = determine_method_and_basis(
        molecule_data,
        functional_override=functional_override,
        basis_set_override=basis_set_override,
    )
    assert result["functional"] == functional_override
    assert result["basis_set"] == basis_set_override


def test_invalid_molecule_data_never_raises():
    for bad_input in [{}, {"valid": False}, None, {"valid": True}]:
        result = determine_method_and_basis(bad_input)
        assert isinstance(result, dict)
