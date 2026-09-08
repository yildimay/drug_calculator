from core.admet import analyze_drug_likeness, print_drug_likeness_summary, LIPINSKI_MAX_TOLERATED_VIOLATIONS

# A large, flexible, polar cyclic peptide-like structure - deliberately
# violates both Lipinski and Veber on multiple counts (MW, HBD, HBA, TPSA,
# rotatable bonds all well over their thresholds).
_LARGE_LIPOPHILIC_SMILES = (
    "CCC1NC(=O)C(C(O)C(C)CC=CC)N(C)C(=O)C(C(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)"
    "C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(C)NC(=O)C(CC(C)C)N(C)C(=O)C(NC(=O)C(CC(C)C)"
    "N(C)C(=O)C(C)NC1=O)C(C)C"
)


def test_invalid_smiles_returns_none():
    assert analyze_drug_likeness("not_a_real_smiles(((") is None


def test_aspirin_passes_both_rules():
    result = analyze_drug_likeness("CC(=O)Oc1ccccc1C(=O)O")

    assert result["lipinski_pass"] is True
    assert result["lipinski_violations"] == []
    assert result["veber_pass"] is True
    assert result["veber_violations"] == []
    assert result["suggestions"] == []
    assert result["molecular_formula"] == "C9H8O4"


def test_large_flexible_molecule_fails_both_rules():
    result = analyze_drug_likeness(_LARGE_LIPOPHILIC_SMILES)

    assert result["molecular_weight"] > 500
    assert result["lipinski_pass"] is False
    assert len(result["lipinski_violations"]) > LIPINSKI_MAX_TOLERATED_VIOLATIONS
    assert result["veber_pass"] is False
    assert len(result["veber_violations"]) > 0
    assert len(result["suggestions"]) == 2  # one for Lipinski, one for Veber


def test_single_lipinski_violation_still_passes():
    # Lipinski's own rule of thumb tolerates exactly one violation. This
    # long-chain ester is highly lipophilic (LogP > 5) but stays within range
    # on molecular weight, H-bond donors, and H-bond acceptors.
    long_chain_ester = "CCOC(=O)CCCCCCCCCCCCCCCCCCCCC"
    result = analyze_drug_likeness(long_chain_ester)

    assert len(result["lipinski_violations"]) == 1
    assert result["lipinski_pass"] is True


def test_summary_does_not_contradict_itself_on_tolerated_violation(caplog):
    # Regression test: a molecule with exactly one Lipinski violation still
    # PASSes (Lipinski's own tolerance), but the summary must say so rather
    # than printing "PASS" immediately followed by an unexplained violation
    # bullet, which reads as self-contradictory.
    import logging
    caplog.set_level(logging.INFO)

    long_chain_ester = "CCOC(=O)CCCCCCCCCCCCCCCCCCCCC"
    result = analyze_drug_likeness(long_chain_ester)
    assert result["lipinski_pass"] is True
    assert len(result["lipinski_violations"]) == 1

    print_drug_likeness_summary(result)

    pass_line = next(r.message for r in caplog.records if "Lipinski's Rule of Five:" in r.message)
    assert "PASS" in pass_line
    assert "tolerated" in pass_line
