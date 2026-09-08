# Examples

Captured input/output pairs for the tool, split into three tiers of
increasing complexity. Each folder has a `command.txt` (what was run) and
an `output.txt` (what it printed) - real output, not written by hand.

## [beginner/](beginner/)
The simplest possible use: one SMILES string, no flags. Shows the full
default report end-to-end on progressively more realistic molecules
(ethanol → benzene → aspirin).

## [intermediate/](intermediate/)
Using the tool's flags: reading a failing drug-likeness screen, overriding
the recommended method, getting full `--verbose` detail, and handling a
real transition-metal complex (carboplatin).

## [expert/](expert/)
Genuinely hard cases: a 113-atom drug molecule (paclitaxel) that stress-
tests the decision tree's branch priority, a heavy-element-containing drug
(levothyroxine) where two branches fire independently, the one example
in this whole set that exercises the "large non-aromatic" DFT branch, and
deliberate error handling for a chemically-ambiguous input.

## Reproducing these yourself

All commands assume you're in the project root with the virtual environment
active and the package installed (`pip install -e .`, see the main
[../README.md](../README.md)):
```bash
cd drug_calculator  # your local clone of this repo
source .venv/bin/activate
drugcalculator --smiles "..."
```

One thing that can make your own output differ slightly from what's
recorded here: the `IUPAC Name:` line depends on a live PubChem lookup. If
you're offline, or PubChem is temporarily unavailable, it'll say
"unavailable" instead. Add `--no-name-lookup` to skip that network call
entirely.

See the main [../README.md](../README.md) for full documentation of every
flag and how the recommendation logic works.
