# Drug Calculator

A CLI tool for early drug-development screening. Give it a molecule as a
SMILES string and it looks up the name, screens basic drug-likeness, and
recommends a DFT functional/basis set to run further calculations with.
It doesn't write calculation input files - it's meant to help you decide
what to run, not to run it for you.

## How it works

One SMILES string in, three things out:

1. **Parsing** (`core/molecule_parser.py`, `core/naming.py`) - RDKit
   validates the SMILES, embeds a 3D conformer, and pulls out formal
   charge, aromaticity, transition metals, heavy elements, etc. Also
   checks for atoms stacked on top of each other, which happens when a
   SMILES uses `.` for fragments that were actually meant to be bonded
   (a metal complex written as loose ligands instead of `Metal(Ligand)`
   bonds is the usual culprit - RDKit won't separate them in space on its
   own). A name lookup goes out to PubChem's API; that's the only network
   call this tool makes, and `--no-name-lookup` turns it off.

2. **Drug-likeness** (`core/admet.py`) - Lipinski's Rule of Five and
   Veber's rule: MW, LogP, H-bond donors/acceptors, polar surface area,
   rotatable bonds. Not a real ADMET or toxicity model, just the
   physicochemical thresholds from those two papers. A fail is worth a
   second look, not a verdict - plenty of real drugs break these rules.

3. **Method selection** (`dft/decision_tree.py`) - a small rule-based
   decision tree over the molecule's own properties, nothing learned or
   fetched from outside:
   - Transition metals → TPSSh + LANL2DZ
   - Aromatic rings → wB97XD
   - Large systems (>50 heavy atoms) → PBE0
   - Otherwise → B3LYP-D3

   Basis set and hardware estimate follow similar rules. Override any of
   it with `--functional`, `--basis-set`, `--charge`, `--multiplicity`.

## Setup

```bash
git clone https://github.com/yildimay/drug_calculator.git
cd drug_calculator
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"   # installs rdkit + pytest/hypothesis, and the drugcalculator command
```

That step also puts a `drugcalculator` command on the venv's PATH, so you
don't need `python main.py` anymore. If you'd rather not install the
package, `pip install -r requirements.txt` plus `python main.py ...`
still works exactly the same.

## Examples

```bash
drugcalculator --smiles "CC(=O)Oc1ccccc1C(=O)O"
drugcalculator --smiles "CC(=O)Oc1ccccc1C(=O)O" --verbose
drugcalculator --smiles "[O-]C(=O)C" --charge -1 --functional M06-2X --basis-set def2-TZVP
drugcalculator --smiles "CCO" --no-name-lookup
```

Aspirin, abbreviated:
```
[1/3] Parsing SMILES string...
      IUPAC Name: 2-acetyloxybenzoic acid (via PubChem)
✓ SMILES parsed successfully

[2/3] Screening drug-likeness...
Lipinski's Rule of Five: PASS
Veber's Rule (oral bioavailability): PASS

[3/3] Determining recommended DFT method...
Functional:           wB97XD
Basis Set:            6-31+G(d,p)
```

See `examples/` for full, real captured input/output across three
difficulty levels.

## Flags

- `--verbose` - full molecular detail (3D geometry, every property) instead of the short summary
- `--no-name-lookup` - skip the PubChem call entirely; useful offline, or if you'd rather not send the SMILES out
- `--charge`, `--multiplicity` - override the formal charge / spin multiplicity RDKit derived from the SMILES
- `--functional`, `--basis-set` - override the auto-picked DFT functional / basis set

Run `drugcalculator --help` for the full list with defaults.

## Tests

```bash
pytest
```

147 tests, four layers:
- unit tests per module
- regression tests against ~17 real drugs (`tests/drug_panel.py`) covering
  every decision-tree/basis-set branch
- property-based tests (`hypothesis`) checking invariants across generated
  inputs rather than fixed examples
- a golden-master snapshot of the whole drug panel's output
  (`tests/golden/drug_panel_snapshot.json`), regenerate deliberately with
  `UPDATE_GOLDEN=1 pytest tests/test_golden_master.py` after reviewing the diff

## Contributing
Fork it, add tests for whatever you change, open a PR.

## License
[MIT](LICENSE)
