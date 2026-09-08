# Expert Examples

Genuinely hard cases: branch-priority edge cases, 100+ atom molecules, and
deliberate error handling.

- **01_large_complex_molecule** - Paclitaxel: 113 atoms, MW comes out to
  853.9 Da (matches the real drug almost exactly). It's both aromatic and
  over the 50-heavy-atom threshold, but gets wB97XD instead of PBE0 since
  the aromatic check runs first in the decision tree. Fails Lipinski hard.
  Hardware estimate scales to 16 processors / 96 GB.
- **02_heavy_element_aromatic** - Levothyroxine, a thyroid hormone with
  four iodine atoms. The functional and basis-set branches fire
  independently here: aromatic rings pick wB97XD, the iodines pick
  LANL2DZ. Also a real example of a single tolerated Lipinski violation
  (MW 776.9 > 500) that still passes overall.
- **03_pure_large_system_branch** - A cyclosporine-like cyclic peptide (85
  heavy atoms, no aromatic rings, no metals). The one case here not
  dominated by the aromatic branch, so it's the cleanest demo of the pure
  PBE0 pick. Long floppy chains usually fail RDKit's 3D embedding; this
  one was picked because it's large and still embeds. Fails both
  drug-likeness rules.
- **04_error_handling_disconnected_fragments** - `N.N.Cl[Pt]Cl` looks like
  it could mean cisplatin, but the `.` separators make it two free NH3
  molecules plus a separate Cl-Pt-Cl fragment, not a bonded complex.
  RDKit doesn't spatially separate disconnected fragments, so this
  silently stacks atoms on top of each other. The tool catches it and
  refuses to proceed instead of reporting a method for a broken geometry.
  Compare with `../intermediate/04_transition_metal_complex`, which is
  properly bonded.
