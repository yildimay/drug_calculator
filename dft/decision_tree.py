"""
Decision tree logic for recommending DFT calculation parameters.

This module analyzes molecular properties (from a parsed SMILES) and recommends
a functional, basis set, and estimated hardware allocation. It does not write
any calculation input files; the output is a recommendation report.
"""

import logging
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


def determine_method_and_basis(
    molecule_data: Dict,
    functional_override: Optional[str] = None,
    basis_set_override: Optional[str] = None,
) -> Dict:
    """
    Determine and recommend an optimal DFT functional, basis set, and estimated hardware resources.

    This function implements a decision tree algorithm that analyzes molecular properties
    to select appropriate calculation parameters.

    Parameters
    ----------
    molecule_data : dict
        Parsed molecule dictionary from parse_smiles() containing:
        - 'valid': Boolean indicating if SMILES is valid
        - 'formal_charge': Formal charge of molecule
        - 'has_transition_metals': Boolean
        - 'has_aromatic_rings': Boolean
        - 'has_heavy_elements': Boolean
        - 'num_atoms_total': Total number of atoms
        - 'num_heavy_atoms': Number of heavy atoms
        - 'atom_list': List of (symbol, x, y, z) tuples
    functional_override : str, optional
        If provided, used instead of the auto-selected functional. The rationale
        will still describe what the heuristic would have picked.
    basis_set_override : str, optional
        If provided, used instead of the auto-selected basis set (and disables
        the auto-detected ECP/diffuse-function flags, since those are tied to
        the specific auto-selected basis set).

    Returns
    -------
    dict
        Dictionary containing:
        - 'functional': DFT functional (TPSSh, PBE0, wB97XD, or B3LYP-D3)
        - 'basis_set': Basis set (LANL2DZ, def2-SVP, 6-31+G(d,p), 6-311G(d,p), or 6-31G(d))
        - 'nprocshared': Number of processors (estimated from atom count)
        - 'mem_per_proc': Memory per processor in GB
        - 'total_mem': Total memory in GB
        - 'dispersion_correction': String indicating if dispersion is included in functional
        - 'ecp_used': Boolean indicating if effective core potential is used
        - 'diffuse_functions': Boolean indicating if diffuse functions are included
        - 'rationale': String explaining the selection decisions
        - 'error': Error message if molecule_data is invalid

    Examples
    --------
    >>> result = determine_method_and_basis(molecule_data)
    >>> print(f"Functional: {result['functional']}, Basis: {result['basis_set']}")
    """
    # Validate input
    if not molecule_data or not molecule_data.get('valid', False):
        return {
            'error': 'Invalid molecule data provided',
            'functional': 'unspecified',
            'basis_set': 'unspecified',
            'nprocshared': 1,
            'mem_per_proc': 2,
            'total_mem': 2
        }

    # Extract molecular properties
    has_transition_metals = molecule_data.get('has_transition_metals', False)
    has_aromatic_rings = molecule_data.get('has_aromatic_rings', False)
    has_heavy_elements = molecule_data.get('has_heavy_elements', False)
    formal_charge = molecule_data.get('formal_charge', 0)
    num_atoms = molecule_data.get('num_atoms_total', 0)
    num_heavy_atoms = molecule_data.get('num_heavy_atoms', 0)
    atom_list = molecule_data.get('atom_list', [])

    # Determine functional (selection and rationale come from the same decision
    # so the two can never drift apart)
    auto_functional, functional_rationale = _select_functional(
        has_transition_metals, has_aromatic_rings, num_heavy_atoms
    )
    functional = functional_override or auto_functional
    if functional_override:
        functional_rationale = (
            f"Functional '{functional_override}' was explicitly requested, overriding the "
            f"auto-selected '{auto_functional}'. Auto-selection reasoning: {functional_rationale}"
        )

    # Determine basis set
    auto_basis_set, ecp_used, diffuse_functions, basis_rationale = _select_basis_set(
        has_heavy_elements,
        formal_charge,
        atom_list,
        num_heavy_atoms
    )
    basis_set = basis_set_override or auto_basis_set
    if basis_set_override:
        basis_rationale = (
            f"Basis set '{basis_set_override}' was explicitly requested, overriding the "
            f"auto-selected '{auto_basis_set}'. Auto-selection reasoning: {basis_rationale}"
        )
        # ECP/diffuse-function flags were derived for the auto-selected basis set;
        # they no longer apply once the user picks their own basis set.
        ecp_used = False
        diffuse_functions = False

    # Estimate hardware resources
    nprocshared, mem_per_proc, total_mem = _estimate_hardware_resources(num_atoms)

    # Generate rationale directly from the reasoning produced during selection
    rationale = f"{functional_rationale} {basis_rationale}"

    # Determine dispersion correction
    dispersion_correction = _get_dispersion_info(functional)

    results = {
        'functional': functional,
        'basis_set': basis_set,
        'nprocshared': nprocshared,
        'mem_per_proc': mem_per_proc,
        'total_mem': total_mem,
        'dispersion_correction': dispersion_correction,
        'ecp_used': ecp_used,
        'diffuse_functions': diffuse_functions,
        'rationale': rationale,
    }

    return results


def _select_functional(
    has_transition_metals: bool, has_aromatic_rings: bool, num_heavy_atoms: int
) -> Tuple[str, str]:
    """
    Select appropriate DFT functional based on molecular composition and size.

    Decision logic:
    1. Transition metals present -> TPSSh (better for metal complexes)
    2. Aromatic rings present -> wB97XD (includes dispersion for pi systems)
    3. Large systems (>50 heavy atoms) -> PBE0 (better scaling for large molecules)
    4. Default -> B3LYP-D3 (balanced for organic molecules)

    Parameters
    ----------
    has_transition_metals : bool
        Whether molecule contains transition metals.
    has_aromatic_rings : bool
        Whether molecule contains aromatic rings.
    num_heavy_atoms : int
        Number of heavy atoms in molecule.

    Returns
    -------
    tuple
        - functional: Selected functional.
        - rationale: Sentence explaining why this functional was chosen.
    """
    if has_transition_metals:
        # TPSSh is generally better for geometry optimization of metal complexes
        return 'TPSSh', (
            "Functional 'TPSSh' selected: Transition metals detected. "
            "This meta-GGA hybrid functional is excellent for metal complexes."
        )
    elif has_aromatic_rings:
        # wB97XD includes dispersion and handles pi-pi interactions well
        return 'wB97XD', (
            "Functional 'wB97XD' selected: Aromatic rings detected. "
            "This range-separated hybrid functional includes dispersion for pi-systems."
        )
    elif num_heavy_atoms > 50:
        # PBE0 scales better for large systems while maintaining accuracy
        return 'PBE0', (
            f"Functional 'PBE0' selected: Large system ({num_heavy_atoms} heavy atoms). "
            "This hybrid functional provides better computational scaling for large molecules."
        )
    else:
        # B3LYP-D3 is the standard workhorse for organic molecules
        return 'B3LYP-D3', (
            "Functional 'B3LYP-D3' selected: Standard organic molecule. "
            "This is the workhorse functional for general chemistry."
        )


def _select_basis_set(
    has_heavy_elements: bool,
    formal_charge: int,
    atom_list: list,
    num_heavy_atoms: int
) -> Tuple[str, bool, bool, str]:
    """
    Select appropriate basis set based on molecular properties.

    Decision logic:
    1. Heavy elements (Period 4+) present -> LANL2DZ or def2-SVP (with ECP)
    2. Charge < 0 or high electronegative atoms -> 6-31+G(d,p) (with diffuse functions)
    3. Small organic molecules (heavy atoms < 20) -> 6-311G(d,p) (high quality)
    4. Larger organic molecules (heavy atoms >= 20) -> 6-31G(d) (balanced)

    Parameters
    ----------
    has_heavy_elements : bool
        Whether molecule contains heavy elements (Period 4+).
    formal_charge : int
        Formal charge of molecule.
    atom_list : list
        List of (symbol, x, y, z) tuples for each atom.
    num_heavy_atoms : int
        Number of heavy atoms in molecule.

    Returns
    -------
    tuple
        - basis_set: Selected basis set string
        - ecp_used: Boolean indicating if ECP is used
        - diffuse_functions: Boolean indicating if diffuse functions are included
        - rationale: Sentence explaining why this basis set was chosen.
    """
    # Check for high electronegative atom count (O, N, F, Cl)
    electronegative_atoms = ['O', 'N', 'F', 'Cl']
    electronegative_count = sum(
        1 for symbol, _, _, _ in atom_list if symbol in electronegative_atoms
    )

    # Rule 1: Heavy elements present
    if has_heavy_elements:
        # LANL2DZ is standard for heavy elements, def2-SVP is more modern alternative
        # Using LANL2DZ for compatibility
        return 'LANL2DZ', True, False, (
            "Basis set 'LANL2DZ' selected: Heavy elements (Period 4+) detected. "
            "Effective core potential basis is necessary for computational efficiency."
        )

    # Rule 2: Anionic species or electron-rich/polarizable system
    if formal_charge < 0 or electronegative_count >= 3:
        # 6-31+G(d,p) includes diffuse functions for anions and electronegative systems.
        # Report whichever condition(s) actually triggered this, rather than always
        # claiming "anionic" even for neutral, electronegative-atom-rich molecules
        # (common for drug-like molecules with several O/N/F/Cl atoms).
        reasons = []
        if formal_charge < 0:
            reasons.append(f"anionic species (charge={formal_charge})")
        if electronegative_count >= 3:
            reasons.append(
                f"{electronegative_count} electronegative atoms (O/N/F/Cl) suggest a "
                "polarizable, lone-pair-rich system"
            )
        reason_text = " and ".join(reasons)
        reason_text = reason_text[0].upper() + reason_text[1:]
        return '6-31+G(d,p)', False, True, (
            f"Basis set '6-31+G(d,p)' selected: {reason_text}. "
            "Diffuse functions added to properly describe the electron density in these cases."
        )

    # Rule 3 & 4: Standard organic molecules - basis size depends on system size
    if num_heavy_atoms < 20:
        # Smaller systems can afford higher quality basis
        return '6-311G(d,p)', False, False, (
            f"Basis set '6-311G(d,p)' selected: Small organic system ({num_heavy_atoms} heavy atoms). "
            "Higher-quality triple-zeta basis suitable for accurate calculations."
        )
    else:
        # Larger systems use more economical basis
        return '6-31G(d)', False, False, (
            f"Basis set '6-31G(d)' selected: Larger system ({num_heavy_atoms} heavy atoms). "
            "Double-zeta basis balances accuracy and computational cost."
        )


def _estimate_hardware_resources(num_atoms: int) -> Tuple[int, float, float]:
    """
    Estimate hardware resources (%nprocshared and %mem) based on atom count.

    General guidelines:
    - Small systems (< 20 atoms): 2-4 processors, 2 GB/core
    - Medium systems (20-50 atoms): 4-8 processors, 2-4 GB/core
    - Large systems (50-100 atoms): 8-16 processors, 4-6 GB/core
    - Very large systems (> 100 atoms): 16+ processors, 6-8 GB/core

    Parameters
    ----------
    num_atoms : int
        Total number of atoms in molecule.

    Returns
    -------
    tuple
        - nprocshared: Number of processors
        - mem_per_proc: Memory per processor in GB
        - total_mem: Total memory in GB
    """
    if num_atoms < 20:
        nprocshared = 2
        mem_per_proc = 2.0
    elif num_atoms < 50:
        nprocshared = 4
        mem_per_proc = 2.5
    elif num_atoms < 100:
        nprocshared = 8
        mem_per_proc = 4.0
    else:
        nprocshared = 16
        mem_per_proc = 6.0

    total_mem = nprocshared * mem_per_proc
    return nprocshared, mem_per_proc, total_mem


def _get_dispersion_info(functional: str) -> str:
    """
    Get information about dispersion correction in functional.

    Parameters
    ----------
    functional : str
        Selected DFT functional.

    Returns
    -------
    str
        Description of dispersion treatment in functional.
    """
    dispersion_map = {
        'TPSSh': 'Built-in dispersion (meta-GGA+HF)',
        'PBE0': 'No explicit dispersion - may need manual D3 correction for accuracy',
        'wB97XD': 'Explicit dispersion (D parametrization)',
        'B3LYP-D3': 'Explicit D3 correction required'
    }
    return dispersion_map.get(functional, 'Standard DFT')


def print_method_summary(method_params: Dict) -> None:
    """
    Print a formatted summary of selected DFT parameters.

    Parameters
    ----------
    method_params : dict
        Dictionary from determine_method_and_basis().
    """
    logger.info("=" * 70)
    logger.info("DFT CALCULATION PARAMETERS")
    logger.info("=" * 70)
    logger.info(f"Functional:           {method_params.get('functional', 'N/A')}")
    logger.info(f"Basis Set:            {method_params.get('basis_set', 'N/A')}")
    logger.info(f"Dispersion:           {method_params.get('dispersion_correction', 'N/A')}")
    logger.info(f"ECP Used:             {method_params.get('ecp_used', False)}")
    logger.info(f"Diffuse Functions:    {method_params.get('diffuse_functions', False)}")
    logger.info("\nHardware Resources:")
    logger.info(f"  Processors:         {method_params.get('nprocshared', 1)}")
    logger.info(f"  Memory/Processor:   {method_params.get('mem_per_proc', 2.0)} GB")
    logger.info(f"  Total Memory:       {method_params.get('total_mem', 2)} GB")
    logger.info("\nRationale:")
    logger.info(f"  {method_params.get('rationale', 'N/A')}")
    logger.info("=" * 70)
