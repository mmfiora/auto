"""
resolve_unusual_aa.py
---------------------
For each peptide in the `peptides` table whose sequence contains 'X':
  1. Fetch the DBAASP API to get the actual unusual AA codes per position.
  2. Store them as JSON in the new `unusual_residues` column
     e.g. {"1": "AIB", "4": "AIB", "8": "AIB"}
  3. Use the PubChem CID (returned by the API) to fetch the full-peptide SMILES.
  4. Compute logP and logD from that SMILES using the existing RDKit pipeline.
  5. UPDATE the `peptides` row with all four fields.

Run from the repo root:
    python dbassp/scripts/resolve_unusual_aa.py
"""

import json
import logging
import os
import sqlite3
import sys
import time

import requests

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.collectors import lipophilicity as lipo

DB_PATH = "dbassp/data/output/dbaasp.sqlite"
DBAASP_URL = "https://dbaasp.org/peptides/{pid}"
PUBCHEM_SMILES_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
    "/property/IsomericSMILES/JSON"
)

DBAASP_DELAY = 0.5   # seconds between DBAASP requests
PUBCHEM_DELAY = 0.3  # seconds between PubChem requests

logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
log = logging.getLogger(__name__)


# ── helpers ───────────────────────────────────────────────────────────────────

def _fetch_json(url: str, timeout: int = 20) -> dict | None:
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
        resp.raise_for_status()
        return resp.json()
    except Exception as exc:
        log.warning(f"Request failed {url}: {exc}")
        return None


def fetch_dbaasp(pid: int) -> dict | None:
    return _fetch_json(DBAASP_URL.format(pid=pid))


def fetch_pubchem_smiles(cid: str) -> str | None:
    # Try fetching both IsomericSMILES and SMILES, fallback if needed
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES,SMILES/JSON"
    data = _fetch_json(url)
    if not data:
        return None
    try:
        props = data["PropertyTable"]["Properties"][0]
        return props.get("IsomericSMILES") or props.get("SMILES")
    except (KeyError, IndexError):
        return None


def get_unusual_residues(api_data: dict) -> dict | None:
    """Return dict mapping position -> AA code, or None if no unusual AAs."""
    entries = api_data.get("unusualAminoAcids") or []
    if not entries:
        return None
    mapping = {}
    for entry in entries:
        pos = entry.get("position")
        code = (entry.get("modificationType") or {}).get("name", "")
        if pos is not None and code:
            mapping[pos] = code
    return mapping if mapping else None


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    # 1. Create junction table if missing
    cur.execute("""
        CREATE TABLE IF NOT EXISTS peptide_unusual_residues (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            peptide_id INTEGER NOT NULL,
            position INTEGER NOT NULL,
            aa_name TEXT NOT NULL,
            FOREIGN KEY (peptide_id) REFERENCES peptides (id),
            FOREIGN KEY (aa_name) REFERENCES unusual_amino_acids (name)
        )
    """)
    conn.commit()

    # 2. Fetch peptides with X in sequence
    cur.execute(
        "SELECT id, sequence, n_terminus, c_terminus FROM peptides WHERE sequence LIKE '%X%'"
    )
    rows = cur.fetchall()
    total = len(rows)
    log.info(f"Found {total} peptides with X in sequence")

    stats = {"unusual_ok": 0, "smiles_ok": 0, "logp_ok": 0, "skipped": 0}

    for i, (pid, seq, n_term, c_term) in enumerate(rows, 1):
        log.info(f"[{i}/{total}] Peptide {pid} — {seq!r}")

        # ── DBAASP API ────────────────────────────────────────────────────────
        api = fetch_dbaasp(pid)
        time.sleep(DBAASP_DELAY)

        if not api:
            log.warning(f"  Skipping {pid}: no API data")
            stats["skipped"] += 1
            continue

        unusual_dict = get_unusual_residues(api)
        if unusual_dict:
            stats["unusual_ok"] += 1
            log.info(f"  Unusual AAs: {unusual_dict}")
        else:
            log.info(f"  No unusual AA data from API")

        # ── PubChem SMILES ─────────────────────────────────────────────────────
        smiles = None
        cid = (api.get("pubChem") or {}).get("cid")
        if cid:
            smiles = fetch_pubchem_smiles(str(cid))
            time.sleep(PUBCHEM_DELAY)
            if smiles:
                stats["smiles_ok"] += 1
                log.info(f"  SMILES from PubChem (CID={cid}): {smiles[:60]}...")
            else:
                log.info(f"  PubChem CID={cid} returned no SMILES")
        else:
            log.info(f"  No PubChem CID")

        # ── logP / logD ────────────────────────────────────────────────────────
        logp = logd = None
        if smiles:
            logp = lipo.calculate_logp(smiles)
            logd = lipo.calculate_logd(
                smiles, sequence=seq, ph=7.0,
                nterminus=n_term, cterminus=c_term
            )
            if logp is not None:
                stats["logp_ok"] += 1
                log.info(f"  logP={logp:.3f}  logD={logd:.3f}")

        # ── UPDATE ─────────────────────────────────────────────────────────────
        cur.execute(
            """UPDATE peptides
               SET smiles           = ?,
                   logp             = ?,
                   logd             = ?
               WHERE id = ?""",
            (smiles, logp, logd, pid),
        )
        
        # Insert unusual residues into junction table
        if unusual_dict:
            # First remove any existing mapping for this peptide to avoid duplicates on re-run
            cur.execute("DELETE FROM peptide_unusual_residues WHERE peptide_id = ?", (pid,))
            for pos, aa_name in unusual_dict.items():
                cur.execute(
                    "INSERT INTO peptide_unusual_residues (peptide_id, position, aa_name) VALUES (?, ?, ?)",
                    (pid, pos, aa_name)
                )

        conn.commit()
        conn.commit()

    conn.close()
    log.info("Done.")
    log.info(f"  unusual_residues populated : {stats['unusual_ok']}/{total}")
    log.info(f"  smiles populated           : {stats['smiles_ok']}/{total}")
    log.info(f"  logP/logD populated        : {stats['logp_ok']}/{total}")
    log.info(f"  skipped (no API)           : {stats['skipped']}")


if __name__ == "__main__":
    main()
