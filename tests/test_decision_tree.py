from dft.decision_tree import determine_method_and_basis


def _base_molecule_data(**overrides):
    data = {
        "valid": True,
        "has_transition_metals": False,
        "has_aromatic_rings": False,
        "has_heavy_elements": False,
        "formal_charge": 0,
        "num_atoms_total": 10,
        "num_heavy_atoms": 5,
        "atom_list": [("C", 0.0, 0.0, 0.0)] * 5,
    }
    data.update(overrides)
    return data


def test_invalid_molecule_data_returns_error():
    result = determine_method_and_basis({"valid": False})
    assert "error" in result


def test_transition_metal_selects_tpssh_and_lanl2dz():
    data = _base_molecule_data(has_transition_metals=True, has_heavy_elements=True)
    result = determine_method_and_basis(data)

    assert result["functional"] == "TPSSh"
    assert result["basis_set"] == "LANL2DZ"
    assert result["ecp_used"] is True


def test_aromatic_selects_wb97xd():
    data = _base_molecule_data(has_aromatic_rings=True)
    result = determine_method_and_basis(data)

    assert result["functional"] == "wB97XD"


def test_large_system_selects_pbe0():
    data = _base_molecule_data(num_heavy_atoms=60)
    result = determine_method_and_basis(data)

    assert result["functional"] == "PBE0"


def test_default_selects_b3lyp_d3():
    data = _base_molecule_data()
    result = determine_method_and_basis(data)

    assert result["functional"] == "B3LYP-D3"


def test_anion_gets_diffuse_functions():
    data = _base_molecule_data(formal_charge=-1)
    result = determine_method_and_basis(data)

    assert result["basis_set"] == "6-31+G(d,p)"
    assert result["diffuse_functions"] is True
    assert "anionic" in result["rationale"].lower()


def test_neutral_electronegative_rich_molecule_not_called_anionic():
    # Regression test: a neutral molecule with >=3 electronegative atoms
    # (e.g. aspirin) also triggers the 6-31+G(d,p)/diffuse-functions rule,
    # but the rationale must not claim it's "anionic" when charge is 0.
    atom_list = [("O", 0.0, 0.0, 0.0)] * 3 + [("C", 0.0, 0.0, 0.0)] * 2
    data = _base_molecule_data(formal_charge=0, atom_list=atom_list)
    result = determine_method_and_basis(data)

    assert result["basis_set"] == "6-31+G(d,p)"
    assert "anionic" not in result["rationale"].lower()
    assert "electronegative" in result["rationale"].lower()


def test_anionic_and_electronegative_rich_mentions_both_reasons():
    atom_list = [("O", 0.0, 0.0, 0.0)] * 3 + [("C", 0.0, 0.0, 0.0)] * 2
    data = _base_molecule_data(formal_charge=-1, atom_list=atom_list)
    result = determine_method_and_basis(data)

    assert "anionic" in result["rationale"].lower()
    assert "electronegative" in result["rationale"].lower()


def test_small_vs_large_organic_basis_set_size():
    small = determine_method_and_basis(_base_molecule_data(num_heavy_atoms=10))
    large = determine_method_and_basis(_base_molecule_data(num_heavy_atoms=30))

    assert small["basis_set"] == "6-311G(d,p)"
    assert large["basis_set"] == "6-31G(d)"


def test_functional_and_basis_overrides_are_applied():
    data = _base_molecule_data(has_aromatic_rings=True)
    result = determine_method_and_basis(
        data, functional_override="M06-2X", basis_set_override="def2-TZVP"
    )

    assert result["functional"] == "M06-2X"
    assert result["basis_set"] == "def2-TZVP"
    # Rationale should mention both the override and what the heuristic would have picked
    assert "M06-2X" in result["rationale"]
    assert "wB97XD" in result["rationale"]


def test_rationale_mentions_selected_functional_and_basis():
    data = _base_molecule_data(has_transition_metals=True, has_heavy_elements=True)
    result = determine_method_and_basis(data)

    assert result["functional"] in result["rationale"]
    assert result["basis_set"] in result["rationale"]


def test_hardware_scales_with_atom_count():
    small = determine_method_and_basis(_base_molecule_data(num_atoms_total=10))
    large = determine_method_and_basis(_base_molecule_data(num_atoms_total=150))

    assert small["nprocshared"] < large["nprocshared"]
    assert small["total_mem"] < large["total_mem"]
