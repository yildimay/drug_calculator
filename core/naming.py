"""
Look up a human-readable IUPAC name for a SMILES string via PubChem.

RDKit has no offline IUPAC name generator, so this queries PubChem's free
public PUG-REST API (no API key required). This is the only network call
in this tool: it sends the molecule's SMILES string to PubChem's servers,
which is worth knowing if you're working with sensitive/proprietary
structures. It fails gracefully (returns None) on any network error, a
non-2xx response (e.g. the structure isn't in PubChem), or a timeout, so
the rest of the tool's output is unaffected if this is unavailable.
"""

import logging
import urllib.error
import urllib.parse
import urllib.request
from typing import Optional

logger = logging.getLogger(__name__)

PUBCHEM_PUG_REST_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
DEFAULT_TIMEOUT_SECONDS = 10.0


def get_iupac_name(smiles: str, timeout: float = DEFAULT_TIMEOUT_SECONDS) -> Optional[str]:
    """
    Fetch the IUPAC name for a SMILES string from PubChem.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule.
    timeout : float, optional
        Network timeout in seconds (default: 10.0).

    Returns
    -------
    str or None
        The IUPAC name if PubChem has one on record, otherwise None (this
        includes the "not found" case, network failures, and timeouts -
        callers should treat None as "unavailable," not as an error).
    """
    encoded_smiles = urllib.parse.quote(smiles, safe="")
    url = f"{PUBCHEM_PUG_REST_BASE}/compound/smiles/{encoded_smiles}/property/IUPACName/TXT"

    request = urllib.request.Request(
        url,
        headers={"User-Agent": "drug-calculator/1.0 (drug-development research tool)"},
    )

    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            name = response.read().decode("utf-8").strip()
            return name or None
    except urllib.error.HTTPError as e:
        if e.code == 404:
            logger.debug(f"PubChem has no record for this structure (404): {smiles}")
        else:
            logger.warning(f"⚠ PubChem lookup failed (HTTP {e.code}): {e.reason}")
        return None
    except urllib.error.URLError as e:
        logger.warning(f"⚠ PubChem lookup unavailable (network error): {e.reason}")
        return None
    except Exception as e:
        logger.warning(f"⚠ PubChem lookup failed unexpectedly: {type(e).__name__}: {e}")
        return None
