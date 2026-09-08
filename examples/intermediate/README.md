# Intermediate Examples

Using the tool's flags: overriding a recommendation, reading a failing
drug-likeness screen, full verbose detail, and a transition-metal complex.

- **01_failing_drug_likeness** - Atorvastatin (Lipitor), a real approved
  drug that fails both Lipinski's Rule of Five and Veber's rule. Shows the
  `Suggestions:` section that only shows up when a rule fails.
- **02_method_override** - `--charge -1 --functional M06-2X --basis-set def2-TZVP`
  on acetate, overriding everything at once. The rationale names both what
  you asked for and what the heuristic would've picked instead.
- **03_verbose_mode** - `--verbose` on aspirin: full 3D detail instead of
  the compact summary.
- **04_transition_metal_complex** - Carboplatin, a real chemotherapy drug,
  properly bonded Pt(II) complex. Gets TPSSh + LANL2DZ.

Try it yourself:
```bash
cd drug_calculator  # your local clone of this repo
drugcalculator --smiles "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O"
```

See the top-level [../README.md](../README.md) for full flag documentation.
