from rdkit import Chem
from rdkit.Chem import AllChem

from core.molecule_parser import parse_smiles, get_molecule_info, _find_atom_clash


def test_parse_valid_smiles_benzene():
    result = parse_smiles("c1ccccc1")

    assert result["valid"] is True
    assert result["formal_charge"] == 0
    assert result["spin_multiplicity"] == 1
    assert result["has_aromatic_rings"] is True
    assert result["has_transition_metals"] is False
    assert result["num_heavy_atoms"] == 6
    assert result["num_atoms_total"] == 12  # 6 C + 6 H
    assert len(result["atom_list"]) == 12


def test_parse_invalid_smiles_returns_error():
    result = parse_smiles("not_a_real_smiles(((")

    assert result["valid"] is False
    assert "error" in result


def test_parse_smiles_detects_transition_metal():
    # Pt(NH3)2Cl2 - cisplatin-like fragment
    result = parse_smiles("[N][Pt]([N])(Cl)Cl")

    assert result["valid"] is True
    assert result["has_transition_metals"] is True
    assert result["has_heavy_elements"] is True


def test_parse_smiles_charge_override():
    baseline = parse_smiles("c1ccccc1")
    overridden = parse_smiles("c1ccccc1", charge_override=-1, multiplicity_override=2)

    assert baseline["formal_charge"] == 0
    assert overridden["formal_charge"] == -1
    assert overridden["spin_multiplicity"] == 2


def test_get_molecule_info_returns_same_as_parse_smiles():
    direct = parse_smiles("CCO")
    via_helper = get_molecule_info("CCO", verbose=False)

    assert direct["formal_charge"] == via_helper["formal_charge"]
    assert direct["num_atoms_total"] == via_helper["num_atoms_total"]


def test_disconnected_fragments_smiles_rejected_as_overlapping():
    # "." separates fragments RDKit's embedder does not spatially separate -
    # this is not real bonded cisplatin (would need N[Pt](N)(Cl)Cl), it's two
    # free NH3 molecules plus a separate Cl-Pt-Cl fragment landing on top of
    # each other.
    result = parse_smiles("N.N.Cl[Pt]Cl")

    assert result["valid"] is False
    assert "overlapping" in result["error"].lower()


def test_ionic_salt_smiles_also_rejected_as_overlapping():
    # Not cisplatin-specific: any multi-fragment SMILES collapses the same way.
    result = parse_smiles("[Na+].[Cl-]")

    assert result["valid"] is False
    assert "overlapping" in result["error"].lower()


def test_normal_single_fragment_molecules_have_no_clash():
    for smiles in ("CCO", "c1ccccc1", "CC(=O)[O-]"):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        assert _find_atom_clash(mol) is None


def test_find_atom_clash_detects_manually_overlapped_atoms():
    mol = Chem.MolFromSmiles("[Na+].[Cl-]")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)

    clash = _find_atom_clash(mol)

    assert clash is not None
    idx_a, idx_b, distance = clash
    assert distance < 0.5
