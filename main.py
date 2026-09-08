"""
Main CLI entry point for the DFT Method Recommender application.

This script provides a command-line interface to:
- Parse a SMILES string and extract molecular features
- Look up a human-readable IUPAC name for the molecule (via PubChem)
- Screen the molecule for basic oral drug-likeness (Lipinski/Veber rules)
- Recommend a DFT functional, basis set, and estimated hardware allocation,
  using a rule-based decision tree over the molecule's own properties

It does not write any calculation input files - the output is a
recommendation report intended to help decide how to run a DFT calculation
(e.g. as part of a drug-development computational pipeline), not a
ready-to-submit Gaussian/ORCA input.
"""

import argparse
import logging
import sys

# Import core functionality
from core.admet import analyze_drug_likeness, print_drug_likeness_summary
from core.molecule_parser import parse_smiles, get_molecule_info
from core.naming import get_iupac_name
from dft.decision_tree import determine_method_and_basis, print_method_summary

logger = logging.getLogger(__name__)


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and return the argument parser for the CLI.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description="""\
DFT Method Recommender - given a SMILES string, looks up its name, screens
its basic drug-likeness, and recommends a DFT functional, basis set, and
estimated hardware allocation for it.

Parses the SMILES with RDKit, embeds and optimizes a 3D geometry, looks up
an IUPAC name via PubChem (network access; skip with --no-name-lookup),
screens the molecule against Lipinski's Rule of Five and Veber's rule (basic
oral drug-likeness heuristics - not a full ADMET/toxicity prediction), then
runs a rule-based decision tree over molecular properties (transition
metals, aromaticity, formal charge, molecule size) to recommend a sensible
DFT functional/basis set. This does not produce a calculation input file -
it's a recommendation report, meant to help decide how to set up a DFT
calculation (e.g. when screening drug candidates) before doing so in your
own Gaussian/ORCA/etc. workflow.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Get a recommended method for benzene
  python main.py --smiles "c1ccccc1"

  # Full molecular detail breakdown
  python main.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --verbose

  # Override the auto-selected functional/basis set and formal charge
  python main.py --smiles "[O-]C(=O)C" --charge -1 --functional M06-2X --basis-set def2-TZVP

  # Skip the PubChem name lookup (fully offline, no network call)
  python main.py --smiles "CCO" --no-name-lookup
        """
    )

    parser.add_argument(
        '--smiles',
        type=str,
        required=True,
        help='SMILES string representing the molecule (e.g., "c1ccccc1" for benzene)'
    )

    parser.add_argument(
        '--charge',
        type=int,
        default=None,
        help='Override the formal charge auto-derived from the SMILES string.'
    )

    parser.add_argument(
        '--multiplicity',
        type=int,
        default=None,
        help='Override the spin multiplicity (2S+1) auto-derived from the SMILES string.'
    )

    parser.add_argument(
        '--functional',
        type=str,
        default=None,
        help='Override the auto-selected DFT functional (e.g., "B3LYP", "M06-2X").'
    )

    parser.add_argument(
        '--basis-set',
        type=str,
        default=None,
        dest='basis_set',
        help='Override the auto-selected basis set (e.g., "6-31G(d)", "def2-TZVP").'
    )

    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Show full molecular details (3D geometry, all properties) '
             'instead of just the recommendation summary.'
    )

    parser.add_argument(
        '--no-name-lookup',
        action='store_true',
        dest='no_name_lookup',
        help='Skip the PubChem IUPAC name lookup (avoids the network call entirely; '
             'useful when offline or when you do not want the SMILES sent externally).'
    )

    return parser


def recommend_dft_method(
    smiles: str,
    verbose: bool = False,
    charge_override: int = None,
    multiplicity_override: int = None,
    functional_override: str = None,
    basis_set_override: str = None,
    skip_name_lookup: bool = False,
) -> bool:
    """
    Parse a SMILES string and report its name, drug-likeness, and a
    recommended DFT functional/basis set.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    verbose : bool, optional
        Show full molecular details (default: False).
    charge_override : int, optional
        Override the formal charge auto-derived from the SMILES string.
    multiplicity_override : int, optional
        Override the spin multiplicity auto-derived from the SMILES string.
    functional_override : str, optional
        Override the auto-selected DFT functional.
    basis_set_override : str, optional
        Override the auto-selected basis set.
    skip_name_lookup : bool, optional
        Skip the PubChem IUPAC name lookup entirely (no network call).

    Returns
    -------
    bool
        True if a recommendation was produced successfully, False otherwise.
    """
    logger.info("\n" + "="*70)
    logger.info("DFT METHOD RECOMMENDER")
    logger.info("="*70)

    # Step 1: Parse SMILES
    logger.info("\n[1/3] Parsing SMILES string...")
    logger.info(f"      SMILES: {smiles}")

    molecule_data = parse_smiles(
        smiles, charge_override=charge_override, multiplicity_override=multiplicity_override
    )

    if not molecule_data or not molecule_data.get('valid', False):
        logger.error("\n✗ Error: Invalid SMILES string or parsing failed.")
        if molecule_data and 'error' in molecule_data:
            logger.error(f"  Details: {molecule_data['error']}")
        return False

    if not skip_name_lookup:
        iupac_name = get_iupac_name(smiles)
        if iupac_name:
            logger.info(f"      IUPAC Name: {iupac_name} (via PubChem)")
        else:
            logger.info("      IUPAC Name: unavailable (not found in PubChem, or lookup failed)")

    if verbose:
        logger.info("\n--- Molecular Information ---")
        get_molecule_info(smiles, verbose=True)
    else:
        logger.info("✓ SMILES parsed successfully")
        logger.info(f"  Molecular Weight: {molecule_data.get('molecular_weight', 0):.2f} g/mol")
        logger.info(f"  Total Atoms: {molecule_data.get('num_atoms_total', 0)}")
        logger.info(f"  Heavy Atoms: {molecule_data.get('num_heavy_atoms', 0)}")
        logger.info(f"  Formal Charge: {molecule_data.get('formal_charge', 0)}")
        logger.info(f"  Spin Multiplicity: {molecule_data.get('spin_multiplicity', 1)}")

    # Step 2: Screen basic drug-likeness (Lipinski/Veber)
    logger.info("\n[2/3] Screening drug-likeness...")

    admet = analyze_drug_likeness(smiles)
    if admet is None:
        logger.warning("⚠ Could not compute drug-likeness properties for this SMILES; skipping.")
    else:
        logger.info("")
        print_drug_likeness_summary(admet)

    # Step 3: Determine recommended DFT parameters
    logger.info("\n[3/3] Determining recommended DFT method...")

    dft_params = determine_method_and_basis(
        molecule_data,
        functional_override=functional_override,
        basis_set_override=basis_set_override,
    )

    if 'error' in dft_params:
        logger.error("\n✗ Error: Could not determine a recommendation.")
        logger.error(f"  Details: {dft_params['error']}")
        return False

    logger.info("")
    print_method_summary(dft_params)
    logger.info("")

    return True


def configure_logging(verbose: bool) -> None:
    """Configure root logging so output looks like plain CLI text (no timestamps/levels)."""
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(message)s',
        stream=sys.stdout,
    )


def main() -> int:
    """
    Main entry point for the CLI application.

    Returns
    -------
    int
        Exit code (0 for success, 1 for failure).
    """
    parser = create_argument_parser()
    args = parser.parse_args()
    configure_logging(args.verbose)

    try:
        success = recommend_dft_method(
            smiles=args.smiles,
            verbose=args.verbose,
            charge_override=args.charge,
            multiplicity_override=args.multiplicity,
            functional_override=args.functional,
            basis_set_override=args.basis_set,
            skip_name_lookup=args.no_name_lookup,
        )
        return 0 if success else 1

    except KeyboardInterrupt:
        logger.warning("\n\n⚠ Process interrupted by user.")
        return 1

    except Exception as e:
        logger.error("\n✗ Unexpected error occurred:")
        logger.error(f"  {type(e).__name__}: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    exit_code = main()
    sys.exit(exit_code)
