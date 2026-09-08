"""
Module for parsing and analyzing SMILES strings using RDKit.

This module provides functionality to validate SMILES strings, generate 3D coordinates,
optimize molecular geometry, and extract molecular properties and features.
"""

import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from typing import Dict, Optional, List, Tuple

logger = logging.getLogger(__name__)


# Transition metals and their atomic numbers
TRANSITION_METALS = {
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30,  # Sc to Zn
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48,  # Y to Cd
    57, 72, 73, 74, 75, 76, 77, 78, 79, 80,  # La, Hf to Hg
    89, 104, 105, 106, 107, 108, 109, 110, 111, 112  # Ac, Rf to Cn
}

# Heavy elements (Period 4 and beyond, excluding transition metals)
HEAVY_ELEMENT_PERIODS = {
    19, 20, 31, 32, 33, 34, 35, 36,  # K, Ca, Ga-Kr (Period 4)
    37, 38, 49, 50, 51, 52, 53, 54,  # Rb, Sr, In-Xe (Period 5)
    55, 56, 81, 82, 83, 84, 85, 86,  # Cs, Ba, Tl-Rn (Period 6)
}


def parse_smiles(
    smiles_str: str,
    charge_override: Optional[int] = None,
    multiplicity_override: Optional[int] = None,
    random_seed: int = 42,
) -> Optional[Dict]:
    """
    Parse and analyze a SMILES string, generating molecular features and 3D coordinates.

    Parameters
    ----------
    smiles_str : str
        SMILES representation of a molecule.
    charge_override : int, optional
        If provided, used instead of the formal charge RDKit derives from the SMILES
        (e.g. for reduced/oxidized species that don't have distinct SMILES notation).
    multiplicity_override : int, optional
        If provided, used instead of the spin multiplicity derived from radical electrons
        (e.g. to request a triplet state on a molecule RDKit parses as a closed-shell singlet).
    random_seed : int, optional
        Random seed for ETKDG 3D embedding (default: 42, for reproducibility).

    Returns
    -------
    dict or None
        A comprehensive dictionary containing:
        - 'smiles': Original SMILES string
        - 'valid': Boolean indicating if SMILES is valid
        - 'formal_charge': Total formal charge of the molecule
        - 'spin_multiplicity': Spin multiplicity (2S+1)
        - 'num_heavy_atoms': Number of heavy atoms (non-hydrogen)
        - 'has_aromatic_rings': Boolean indicating presence of aromatic rings
        - 'has_transition_metals': Boolean indicating presence of transition metals
        - 'has_heavy_elements': Boolean indicating presence of heavy elements (Period 4+)
        - 'molecular_weight': Molecular weight in g/mol
        - 'num_atoms_total': Total number of atoms (including hydrogens)
        - 'num_heavy_atoms_with_h': Number of non-H atoms in structured form
        - 'xyz_coordinates': XYZ coordinate string in standard format
        - 'atom_list': List of (atom_symbol, x, y, z) tuples

        Returns None if SMILES validation fails or 3D coordinate generation fails.

    Notes
    -----
    - Explicit hydrogens are added to the molecule
    - 3D coordinates are generated using ETKDG algorithm
    - Molecular geometry is optimized using MMFF94 force field
    """
    # Validate SMILES
    mol = Chem.MolFromSmiles(smiles_str)
    if mol is None:
        return {
            'smiles': smiles_str,
            'valid': False,
            'error': 'Invalid SMILES string'
        }

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates using ETKDG
    try:
        embed_result = AllChem.EmbedMolecule(mol, randomSeed=random_seed)
    except Exception as e:
        return {
            'smiles': smiles_str,
            'valid': False,
            'error': f'Failed to generate 3D coordinates: {str(e)}'
        }

    if embed_result == -1:
        return {
            'smiles': smiles_str,
            'valid': False,
            'error': 'ETKDG embedding failed to converge on a 3D conformer for this molecule'
        }

    # Optimize geometry using MMFF94
    try:
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is not None:
            AllChem.MMFFOptimizeMolecule(mol, props)
    except Exception as e:
        # If MMFF94 fails, try UFF as fallback
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as uff_e:
            return {
                'smiles': smiles_str,
                'valid': False,
                'error': f'Failed to optimize geometry: {str(e)}, UFF fallback also failed: {str(uff_e)}'
            }

    # RDKit's embedder does not translate disconnected fragments apart from
    # each other - a SMILES with "." (salts, hydrates, or a metal complex
    # mistakenly written as unbonded ligands instead of Metal(Ligand) bonds)
    # can come out with two fragments occupying the exact same coordinates.
    # Catch that here rather than silently writing a degenerate geometry.
    clash = _find_atom_clash(mol)
    if clash is not None:
        idx_a, idx_b, distance = clash
        return {
            'smiles': smiles_str,
            'valid': False,
            'error': (
                f"Embedded geometry has overlapping/near-coincident atoms "
                f"(atom {idx_a} and atom {idx_b}, {distance:.3f} Angstrom apart). "
                f"This usually means the SMILES uses \".\" to separate fragments "
                f"that were actually meant to be bonded (e.g. a metal complex "
                f"written as unbonded ligands like \"N.N.Cl[Pt]Cl\" instead of "
                f"\"N[Pt](N)(Cl)Cl\"), since RDKit's 3D embedder does not "
                f"automatically separate disconnected fragments in space. "
                f"Check the SMILES connectivity."
            )
        }

    # Extract molecular features
    formal_charge = charge_override if charge_override is not None else Chem.GetFormalCharge(mol)
    spin_multiplicity = (
        multiplicity_override if multiplicity_override is not None else _calculate_spin_multiplicity(mol)
    )
    num_heavy_atoms = Descriptors.HeavyAtomCount(mol)
    has_aromatic_rings = _has_aromatic_rings(mol)
    has_transition_metals = _has_transition_metals(mol)
    has_heavy_elements = _has_heavy_elements(mol)

    # Extract molecular weight
    molecular_weight = Descriptors.MolWt(mol)

    # Extract XYZ coordinates
    xyz_string, atom_list = _extract_xyz_coordinates(mol)
    num_atoms_total = mol.GetNumAtoms()

    # Build results dictionary
    results = {
        'smiles': smiles_str,
        'valid': True,
        'formal_charge': formal_charge,
        'spin_multiplicity': spin_multiplicity,
        'num_heavy_atoms': num_heavy_atoms,
        'has_aromatic_rings': has_aromatic_rings,
        'has_transition_metals': has_transition_metals,
        'has_heavy_elements': has_heavy_elements,
        'molecular_weight': molecular_weight,
        'num_atoms_total': num_atoms_total,
        'num_heavy_atoms_with_h': num_heavy_atoms,
        'xyz_coordinates': xyz_string,
        'atom_list': atom_list
    }

    return results


# No real covalent bond is this short (H-H in H2, the shortest common one, is
# ~0.74 Angstrom); anything under this is a degenerate/overlapping embedding,
# not just a strained geometry.
ATOM_CLASH_THRESHOLD_ANGSTROM = 0.5


def _find_atom_clash(
    mol: Chem.Mol, threshold: float = ATOM_CLASH_THRESHOLD_ANGSTROM
) -> Optional[Tuple[int, int, float]]:
    """
    Check an embedded conformer for near-coincident (overlapping) atoms.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule with a 3D conformer.
    threshold : float, optional
        Minimum acceptable distance in Angstrom between any two atoms
        (default: ATOM_CLASH_THRESHOLD_ANGSTROM).

    Returns
    -------
    tuple or None
        (atom_idx_a, atom_idx_b, distance) for the first clashing pair found,
        or None if every pair of atoms is farther apart than the threshold.
    """
    conf = mol.GetConformer()
    positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]

    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            distance = positions[i].Distance(positions[j])
            if distance < threshold:
                return i, j, distance

    return None


def _calculate_spin_multiplicity(mol: Chem.Mol) -> int:
    """
    Calculate spin multiplicity (2S+1) from the number of radical electrons.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object.

    Returns
    -------
    int
        Spin multiplicity (2S+1).
    """
    num_radicals = 0
    for atom in mol.GetAtoms():
        num_radicals += atom.GetNumRadicalElectrons()

    # For a molecule with n unpaired electrons, spin multiplicity = 2S + 1 = n + 1
    spin_multiplicity = num_radicals + 1
    return spin_multiplicity


def _has_aromatic_rings(mol: Chem.Mol) -> bool:
    """
    Check if molecule contains aromatic rings.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object.

    Returns
    -------
    bool
        True if molecule contains aromatic atoms, False otherwise.
    """
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            return True
    return False


def _has_transition_metals(mol: Chem.Mol) -> bool:
    """
    Check if molecule contains transition metals.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object.

    Returns
    -------
    bool
        True if molecule contains transition metal atoms, False otherwise.
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in TRANSITION_METALS:
            return True
    return False


def _has_heavy_elements(mol: Chem.Mol) -> bool:
    """
    Check if molecule contains heavy elements from Period 4 and beyond.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object.

    Returns
    -------
    bool
        True if molecule contains heavy elements, False otherwise.
    """
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        # Check if it's a heavy element or transition metal beyond Period 3
        if atomic_num in HEAVY_ELEMENT_PERIODS or (atomic_num in TRANSITION_METALS):
            return True
    return False


def _extract_xyz_coordinates(mol: Chem.Mol) -> Tuple[str, List[Tuple[str, float, float, float]]]:
    """
    Extract XYZ coordinates from molecule in standard format.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule object with 3D coordinates.

    Returns
    -------
    tuple
        - xyz_string: Multi-line string with XYZ coordinates in standard format
        - atom_list: List of (symbol, x, y, z) tuples for each atom
    """
    conf = mol.GetConformer()
    lines = []
    atom_list = []

    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        symbol = atom.GetSymbol()
        x, y, z = pos.x, pos.y, pos.z

        # Format position with 6 decimal places, common in XYZ files
        line = f"{symbol:2s}    {x:12.6f}    {y:12.6f}    {z:12.6f}"
        lines.append(line)
        atom_list.append((symbol, x, y, z))

    xyz_string = "\n".join(lines)
    return xyz_string, atom_list


def get_molecule_info(smiles_str: str, verbose: bool = False) -> Optional[Dict]:
    """
    Convenience function to get molecule information with optional verbose output.

    Parameters
    ----------
    smiles_str : str
        SMILES representation of a molecule.
    verbose : bool, optional
        If True, print detailed information about the molecule (default: False).

    Returns
    -------
    dict or None
        Same as parse_smiles() return value.
    """
    results = parse_smiles(smiles_str)

    if verbose and results:
        if results.get('valid', False):
            logger.info(f"SMILES: {results['smiles']}")
            logger.info(f"Molecular Weight: {results['molecular_weight']:.2f} g/mol")
            logger.info(f"Formal Charge: {results['formal_charge']}")
            logger.info(f"Spin Multiplicity: {results['spin_multiplicity']}")
            logger.info(f"Heavy Atoms: {results['num_heavy_atoms']}")
            logger.info(f"Total Atoms: {results['num_atoms_total']}")
            logger.info(f"Aromatic Rings: {results['has_aromatic_rings']}")
            logger.info(f"Transition Metals: {results['has_transition_metals']}")
            logger.info(f"Heavy Elements (Period 4+): {results['has_heavy_elements']}")
            logger.info("\nXYZ Coordinates:")
            logger.info(results['xyz_coordinates'])
        else:
            logger.error(f"Error: {results.get('error', 'Unknown error')}")

    return results
