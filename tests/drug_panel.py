"""
A curated panel of real, verified drug/reference-compound SMILES, shared
across regression, property-based, and golden-master tests.

Every entry's expected values were computed by actually running the entry
through core/molecule_parser.py, core/admet.py, and dft/decision_tree.py and
reading back the result (not hand-calculated/guessed) - see the docstring
of each test module that consumes this panel for how it's used. This is a
hand-curated panel rather than a full external dataset (e.g. the original
1997 Lipinski Rule-of-Five WDI set of 2245 compounds): this project has no
network access to fetch such a dataset, and a project this size doesn't
need thousands of compounds for good branch coverage - it needs each
decision-tree/ADME branch exercised by at least one real molecule, which
this panel does.

Deliberately covers, across the panel:
- All four DFT functional branches (transition metal / aromatic / large
  non-aromatic / default)
- All four basis-set branches (heavy element / anion-or-electronegative /
  small organic / large organic)
- Lipinski/Veber PASS and FAIL cases
- Neutral, anionic, and zwitterionic charge states
- A transition metal (Pt) and a heavy main-group element (I)
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class DrugPanelEntry:
    name: str
    smiles: str
    expected_functional: str
    expected_basis_set: str
    expected_has_transition_metals: bool
    expected_has_aromatic_rings: bool
    expected_has_heavy_elements: bool
    expected_formal_charge: int
    expected_lipinski_pass: bool
    expected_veber_pass: bool
    notes: str = ""


DRUG_PANEL = [
    DrugPanelEntry(
        name="aspirin",
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        expected_functional="wB97XD",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
    ),
    DrugPanelEntry(
        name="ibuprofen",
        smiles="CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        expected_functional="wB97XD",
        expected_basis_set="6-311G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Small organic (heavy<20) -> triple-zeta basis branch, distinct from aspirin's diffuse-function branch",
    ),
    DrugPanelEntry(
        name="caffeine",
        smiles="Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        expected_functional="wB97XD",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
    ),
    DrugPanelEntry(
        name="paracetamol",
        smiles="CC(=O)Nc1ccc(O)cc1",
        expected_functional="wB97XD",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
    ),
    DrugPanelEntry(
        name="warfarin",
        smiles="CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
        expected_functional="wB97XD",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
    ),
    DrugPanelEntry(
        name="metformin",
        smiles="CN(C)C(=N)NC(=N)N",
        expected_functional="B3LYP-D3",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=False,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Non-aromatic default functional branch; 5 N atoms trigger the electronegative-rich basis rule despite neutral charge",
    ),
    DrugPanelEntry(
        name="glycine_zwitterion",
        smiles="[NH3+]CC(=O)[O-]",
        expected_functional="B3LYP-D3",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=False,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Net formal charge 0 despite an internal +1/-1 pair - basis rule fires on electronegative atom count (N,O,O), not charge",
    ),
    DrugPanelEntry(
        name="acetate_anion",
        smiles="CC(=O)[O-]",
        expected_functional="B3LYP-D3",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=False,
        expected_has_heavy_elements=False,
        expected_formal_charge=-1,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Genuinely anionic (charge=-1), unlike the zwitterion above",
    ),
    DrugPanelEntry(
        name="imatinib",
        smiles="Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
        expected_functional="wB97XD",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Realistic kinase-inhibitor complexity (37 heavy atoms) that still satisfies Ro5 - contrast case against the FAIL entries below",
    ),
    DrugPanelEntry(
        name="atorvastatin",
        smiles="CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
        expected_functional="wB97XD",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=False,
        expected_veber_pass=False,
        notes="Real Ro5-breaking oral drug: MW 558.7 + LogP 6.31 (Lipinski), 12 rotatable bonds (Veber)",
    ),
    DrugPanelEntry(
        name="paclitaxel",
        smiles=(
            "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(c1ccccc1)"
            "NC(=O)c1ccccc1)O)O)OC(=O)c1ccccc1)(CO4)OC(=O)C)O)C)OC(=O)C"
        ),
        expected_functional="wB97XD",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=False,
        expected_veber_pass=False,
        notes=(
            "62 heavy atoms AND aromatic - demonstrates the aromatic branch winning over the "
            "large-system branch, since functional selection checks aromaticity before size. "
            "MW comes out to 853.9 Da, matching paclitaxel's real molecular weight."
        ),
    ),
    DrugPanelEntry(
        name="carboplatin",
        smiles="C1CC2(C1)C(=O)O[Pt](N)(N)OC2=O",
        expected_functional="TPSSh",
        expected_basis_set="LANL2DZ",
        expected_has_transition_metals=True,
        expected_has_aromatic_rings=False,
        expected_has_heavy_elements=True,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Real bonded Pt(II) chemotherapy drug, not the disconnected-fragment form",
    ),
    DrugPanelEntry(
        name="cisplatin_bonded",
        smiles="N[Pt](N)(Cl)Cl",
        expected_functional="TPSSh",
        expected_basis_set="LANL2DZ",
        expected_has_transition_metals=True,
        expected_has_aromatic_rings=False,
        expected_has_heavy_elements=True,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Properly bonded (contrast against the disconnected 'N.N.Cl[Pt]Cl' clash-rejection test case)",
    ),
    DrugPanelEntry(
        name="levothyroxine",
        smiles="Ic1cc(CC(N)C(=O)O)cc(I)c1Oc1cc(I)c(O)c(I)c1",
        expected_functional="wB97XD",
        expected_basis_set="LANL2DZ",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=True,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes=(
            "Aromatic (functional branch) AND heavy-element iodine (basis branch) firing "
            "independently. MW 776.9 Da > 500 is a single tolerated Lipinski violation - still PASS."
        ),
    ),
    DrugPanelEntry(
        name="testosterone",
        smiles="CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C",
        expected_functional="B3LYP-D3",
        expected_basis_set="6-31G(d)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=False,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Non-aromatic steroid nucleus (fused saturated rings, not aromatic); 21 heavy atoms -> 'larger organic' basis branch",
    ),
    DrugPanelEntry(
        name="penicillin_g",
        smiles="CC1(C(N2C(S1)C(C2=O)NC(=O)Cc3ccccc3)C(=O)O)C",
        expected_functional="wB97XD",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=True,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=True,
        expected_veber_pass=True,
        notes="Contains sulfur (period 3) - confirms S does not trigger the heavy-element (period 4+) rule",
    ),
    DrugPanelEntry(
        name="cyclosporine_like",
        smiles=(
            "CCC1NC(=O)C(C(O)C(C)CC=CC)N(C)C(=O)C(C(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)"
            "C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(C)NC(=O)C(CC(C)C)N(C)C(=O)C(NC(=O)C(CC(C)C)"
            "N(C)C(=O)C(C)NC1=O)C(C)C"
        ),
        expected_functional="PBE0",
        expected_basis_set="6-31+G(d,p)",
        expected_has_transition_metals=False,
        expected_has_aromatic_rings=False,
        expected_has_heavy_elements=False,
        expected_formal_charge=0,
        expected_lipinski_pass=False,
        expected_veber_pass=False,
        notes=(
            "The only panel entry that exercises the pure 'large (>50 heavy atoms) AND "
            "non-aromatic' PBE0 branch - large floppy chains often fail ETKDG embedding, but "
            "this cyclic peptide-like structure embeds reliably."
        ),
    ),
]
