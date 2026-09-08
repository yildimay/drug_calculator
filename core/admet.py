"""
Rule-based drug-likeness screening: Lipinski's Rule of Five and Veber's rule.

Scope note: these are physicochemical heuristics correlated with oral
absorption/bioavailability, not a real ADMET (Absorption, Distribution,
Metabolism, Excretion, Toxicity) prediction - no metabolic pathway,
clearance rate, or toxicity endpoint is modeled, and there's no structural-
alert/toxicophore screening. A "pass" means nothing obviously wrong with
basic oral drug-likeness, not a safety or efficacy guarantee. A "fail" is
worth a second look, not a verdict - plenty of approved drugs (CNS drugs,
natural products, biologics-adjacent compounds especially) break these
rules on purpose.
"""

import logging
from typing import Dict, List, Optional

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

logger = logging.getLogger(__name__)

# Lipinski's own rule of thumb: absorption/permeability problems are more
# likely when *more than one* of the four criteria is violated - a single
# violation is common among approved oral drugs and isn't itself a red flag.
LIPINSKI_MAX_TOLERATED_VIOLATIONS = 1


def analyze_drug_likeness(smiles: str) -> Optional[Dict]:
    """
    Screen a molecule against Lipinski's Rule of Five and Veber's rule.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    dict or None
        None if the SMILES is invalid. Otherwise a dict containing:
        - 'molecular_weight', 'logp', 'num_h_donors', 'num_h_acceptors',
          'tpsa', 'num_rotatable_bonds': raw descriptor values
        - 'molecular_formula': Hill-notation molecular formula
        - 'lipinski_violations': list of human-readable violation strings
        - 'lipinski_pass': bool (True if violations <= 1)
        - 'veber_violations': list of human-readable violation strings
        - 'veber_pass': bool (True if no Veber violations)
        - 'suggestions': list of human-readable suggestion strings (empty if
          both rules pass cleanly)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    molecular_weight = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    num_h_donors = Lipinski.NumHDonors(mol)
    num_h_acceptors = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    num_rotatable_bonds = Lipinski.NumRotatableBonds(mol)
    molecular_formula = rdMolDescriptors.CalcMolFormula(mol)

    lipinski_violations = []
    if molecular_weight > 500:
        lipinski_violations.append(f"Molecular weight {molecular_weight:.1f} Da > 500 Da")
    if logp > 5:
        lipinski_violations.append(f"LogP {logp:.2f} > 5")
    if num_h_donors > 5:
        lipinski_violations.append(f"H-bond donors {num_h_donors} > 5")
    if num_h_acceptors > 10:
        lipinski_violations.append(f"H-bond acceptors {num_h_acceptors} > 10")
    lipinski_pass = len(lipinski_violations) <= LIPINSKI_MAX_TOLERATED_VIOLATIONS

    veber_violations = []
    if tpsa > 140:
        veber_violations.append(f"TPSA {tpsa:.1f} Å² > 140 Å²")
    if num_rotatable_bonds > 10:
        veber_violations.append(f"Rotatable bonds {num_rotatable_bonds} > 10")
    veber_pass = len(veber_violations) == 0

    suggestions: List[str] = []
    if not lipinski_pass:
        suggestions.append(
            "Lipinski Rule of Five: more than one violation present "
            f"({', '.join(lipinski_violations)}) - oral absorption/permeability "
            "may be reduced. Consider reducing molecular weight and/or "
            "lipophilicity, or adding polar groups, if oral bioavailability is "
            "a design goal. (A single violation alone is common among approved "
            "oral drugs and isn't a strong signal on its own.)"
        )
    if not veber_pass:
        suggestions.append(
            "Veber's rule violated "
            f"({', '.join(veber_violations)}) - high polar surface area and/or "
            "many rotatable bonds correlate with poor oral bioavailability in "
            "rat models. Consider reducing flexible-chain length or polar "
            "surface area if oral dosing is the goal."
        )

    return {
        'molecular_weight': molecular_weight,
        'logp': logp,
        'num_h_donors': num_h_donors,
        'num_h_acceptors': num_h_acceptors,
        'tpsa': tpsa,
        'num_rotatable_bonds': num_rotatable_bonds,
        'molecular_formula': molecular_formula,
        'lipinski_violations': lipinski_violations,
        'lipinski_pass': lipinski_pass,
        'veber_violations': veber_violations,
        'veber_pass': veber_pass,
        'suggestions': suggestions,
    }


def print_drug_likeness_summary(admet: Dict) -> None:
    """
    Print a formatted drug-likeness summary.

    Parameters
    ----------
    admet : dict
        Dictionary from analyze_drug_likeness().
    """
    logger.info("=" * 70)
    logger.info("DRUG-LIKENESS SCREENING (not a full ADMET/toxicity prediction)")
    logger.info("=" * 70)
    logger.info(f"Molecular Formula:    {admet['molecular_formula']}")
    logger.info(f"Molecular Weight:     {admet['molecular_weight']:.2f} Da")
    logger.info(f"LogP:                 {admet['logp']:.2f}")
    logger.info(f"H-Bond Donors:        {admet['num_h_donors']}")
    logger.info(f"H-Bond Acceptors:     {admet['num_h_acceptors']}")
    logger.info(f"TPSA:                 {admet['tpsa']:.2f} Å²")
    logger.info(f"Rotatable Bonds:      {admet['num_rotatable_bonds']}")

    lipinski_status = "PASS" if admet['lipinski_pass'] else "FAIL"
    if admet['lipinski_pass'] and admet['lipinski_violations']:
        # Lipinski's own rule of thumb tolerates a single violation - make
        # that explicit here, otherwise "PASS" immediately followed by a
        # listed violation reads as a contradiction.
        n = len(admet['lipinski_violations'])
        lipinski_status += f" ({n} violation{'s' if n != 1 else ''} tolerated)"
    logger.info(f"\nLipinski's Rule of Five: {lipinski_status}")
    if admet['lipinski_violations']:
        for violation in admet['lipinski_violations']:
            logger.info(f"  - {violation}")

    logger.info("\nVeber's Rule (oral bioavailability): " + ("PASS" if admet['veber_pass'] else "FAIL"))
    if admet['veber_violations']:
        for violation in admet['veber_violations']:
            logger.info(f"  - {violation}")

    if admet['suggestions']:
        logger.info("\nSuggestions:")
        for suggestion in admet['suggestions']:
            logger.info(f"  - {suggestion}")

    logger.info("=" * 70)
