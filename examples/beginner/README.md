# Beginner Examples

The simplest possible use: one SMILES string, no flags. Each folder has a
`command.txt` and the real `output.txt` it produced.

- **01_ethanol** - smallest possible input, the default (non-aromatic,
  small) branch: B3LYP-D3 / 6-311G(d,p).
- **02_benzene** - simplest aromatic case, shows the aromatic-ring branch: wB97XD.
- **03_aspirin** - a real drug you'd recognize, PubChem name lookup, clean
  drug-likeness PASS.

Try it yourself:
```bash
cd drug_calculator  # your local clone of this repo
drugcalculator --smiles "CCO"
```

Note: the "IUPAC Name" line needs internet access (queries PubChem). Offline,
it'll say "unavailable" instead - everything else still works. Add
`--no-name-lookup` to skip the call entirely.
