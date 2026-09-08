"""
Golden-master (snapshot) tests: the full recommendation report for a fixed
panel of real drugs is computed and diffed against a frozen expected
snapshot (tests/golden/drug_panel_snapshot.json). This catches unintended
behavior drift that unit tests focused on one function at a time can miss -
e.g. a refactor that changes how two correct-in-isolation pieces combine.

Deliberately excluded from the snapshot: raw 3D coordinates (xyz_coordinates/
atom_list) and free-text rationale wording. Embedded geometry can legitimately
shift between RDKit versions (different ETKDG/MMFF implementations), and
rationale is prose that can be reworded without being wrong - neither is
what this test is meant to guard. Floating-point descriptor values are
rounded before comparison for the same reason (avoid failing on harmless
last-digit differences across platforms/library versions) while still
catching a real behavior change.

To intentionally update the snapshot after a verified, correct change:
    UPDATE_GOLDEN=1 pytest tests/test_golden_master.py
then review the resulting diff in tests/golden/drug_panel_snapshot.json
before committing it - a passing "update" run is not itself a correctness
check, only a way to record new expected values.
"""

import json
import os

import pytest

from core.admet import analyze_drug_likeness
from core.molecule_parser import parse_smiles
from dft.decision_tree import determine_method_and_basis
from tests.drug_panel import DRUG_PANEL

GOLDEN_FILE = os.path.join(os.path.dirname(__file__), "golden", "drug_panel_snapshot.json")
UPDATE_GOLDEN = os.environ.get("UPDATE_GOLDEN") == "1"


def _compute_snapshot(smiles: str) -> dict:
    """Compute the stable (version/wording-independent) parts of the full report for one SMILES."""
    molecule_data = parse_smiles(smiles)
    dft_params = determine_method_and_basis(molecule_data)
    admet = analyze_drug_likeness(smiles)

    return {
        "molecule": {
            "formal_charge": molecule_data["formal_charge"],
            "spin_multiplicity": molecule_data["spin_multiplicity"],
            "num_heavy_atoms": molecule_data["num_heavy_atoms"],
            "num_atoms_total": molecule_data["num_atoms_total"],
            "has_aromatic_rings": molecule_data["has_aromatic_rings"],
            "has_transition_metals": molecule_data["has_transition_metals"],
            "has_heavy_elements": molecule_data["has_heavy_elements"],
            "molecular_weight": round(molecule_data["molecular_weight"], 2),
        },
        "dft": {
            "functional": dft_params["functional"],
            "basis_set": dft_params["basis_set"],
            "ecp_used": dft_params["ecp_used"],
            "diffuse_functions": dft_params["diffuse_functions"],
            "dispersion_correction": dft_params["dispersion_correction"],
            "nprocshared": dft_params["nprocshared"],
            "mem_per_proc": dft_params["mem_per_proc"],
            "total_mem": dft_params["total_mem"],
        },
        "admet": {
            "molecular_formula": admet["molecular_formula"],
            "molecular_weight": round(admet["molecular_weight"], 2),
            "logp": round(admet["logp"], 2),
            "num_h_donors": admet["num_h_donors"],
            "num_h_acceptors": admet["num_h_acceptors"],
            "tpsa": round(admet["tpsa"], 2),
            "num_rotatable_bonds": admet["num_rotatable_bonds"],
            "lipinski_pass": admet["lipinski_pass"],
            "veber_pass": admet["veber_pass"],
            "num_lipinski_violations": len(admet["lipinski_violations"]),
            "num_veber_violations": len(admet["veber_violations"]),
        },
    }


def _load_golden() -> dict:
    if not os.path.isfile(GOLDEN_FILE):
        return {}
    with open(GOLDEN_FILE, "r", encoding="utf-8") as f:
        return json.load(f)


def _write_golden(snapshot: dict) -> None:
    os.makedirs(os.path.dirname(GOLDEN_FILE), exist_ok=True)
    with open(GOLDEN_FILE, "w", encoding="utf-8") as f:
        json.dump(snapshot, f, indent=2, sort_keys=True)
        f.write("\n")


@pytest.fixture(scope="module")
def current_snapshot():
    return {entry.name: _compute_snapshot(entry.smiles) for entry in DRUG_PANEL}


def test_update_golden_master_if_requested(current_snapshot):
    if not UPDATE_GOLDEN:
        pytest.skip("Set UPDATE_GOLDEN=1 to (re)generate the golden snapshot file.")
    _write_golden(current_snapshot)


@pytest.mark.parametrize("entry", DRUG_PANEL, ids=[e.name for e in DRUG_PANEL])
def test_matches_golden_snapshot(entry, current_snapshot):
    if UPDATE_GOLDEN:
        pytest.skip("Running in UPDATE_GOLDEN mode; snapshot comparison skipped.")

    golden = _load_golden()
    assert entry.name in golden, (
        f"No golden snapshot recorded for '{entry.name}' - run with UPDATE_GOLDEN=1 to create one."
    )
    assert current_snapshot[entry.name] == golden[entry.name]
