"""
scrape_unusual_aa.py
--------------------
Scrapes the full list of unusual amino acids from the DBAASP REST API
and populates the `unusual_amino_acids` table in the SQLite database.

API endpoint (discovered from DBAASP JS source):
    GET https://dbaasp.org/lookups/unusualAminoAcid?limit=2000

Each record:
    {
        "name":        "Orn",             # DBAASP abbreviation / code
        "description": "Ornithine; C5H12N2O2"  # full name + formula
    }

The description field uses two formats:
    "Full Name; MolFormula"
    "Full Name; Synonym; MolFormula"

The script:
    1. Fetches data from the API.
    2. Parses name, full_name, synonyms, and molecular_formula from description.
    3. Creates the `unusual_amino_acids` table if it does not exist.
    4. Upserts all records (INSERT OR REPLACE) — safe to re-run.
    5. Prints a summary.

Usage:
    cd dbassp
    python scripts/scrape_unusual_aa.py [--db PATH]
"""

import argparse
import logging
import re
import sqlite3
import sys
import requests

logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
log = logging.getLogger("scrape_unusual_aa")

API_URL  = "https://dbaasp.org/lookups/unusualAminoAcid"
API_PARAMS = {"limit": 2000}
API_HEADERS = {
    "Accept": "application/json",
    "User-Agent": "Mozilla/5.0",
}
DEFAULT_DB = "data/output/dbaasp.sqlite"


# ── Schema ────────────────────────────────────────────────────────────────────

CREATE_TABLE_SQL = """
CREATE TABLE IF NOT EXISTS unusual_amino_acids (
    name             TEXT PRIMARY KEY,   -- DBAASP abbreviation/code
    full_name        TEXT,               -- primary full name
    synonyms         TEXT,               -- semicolon-separated synonyms (may be empty)
    molecular_formula TEXT,              -- e.g. C5H12N2O2
    description_raw  TEXT               -- full original description string
)
"""

UPSERT_SQL = """
INSERT OR REPLACE INTO unusual_amino_acids
    (name, full_name, synonyms, molecular_formula, description_raw)
VALUES (?, ?, ?, ?, ?)
"""


# ── Parsing ───────────────────────────────────────────────────────────────────

_FORMULA_RE = re.compile(r"\b[A-Z][a-z]?(?:\d+)?(?:[A-Z][a-z]?(?:\d+)?)+\d*\b")


def _looks_like_formula(token: str) -> bool:
    """Heuristic: a molecular formula contains C and either H or N."""
    t = token.strip()
    return bool(re.match(r"^[A-Z][A-Za-z0-9]*$", t)) and "C" in t and any(x in t for x in ("H", "N", "O"))


def _split_outside_parens(text: str, sep: str = ",") -> list[str]:
    """Split `text` on `sep` only when not inside parentheses."""
    parts, current, depth = [], [], 0
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(ch)
    parts.append("".join(current).strip())
    return [p for p in parts if p]


def parse_description(raw: str) -> tuple[str, str, str]:
    """
    Parse the DBAASP description string.

    Returns (full_name, synonyms_str, molecular_formula).

    Examples:
        "Ornithine; C5H12N2O2"
            → ("Ornithine", "", "C5H12N2O2")
        "(1R,2R)-2-aminocyclopentane carboxylic acid, C6H11NO2"
            → ("(1R,2R)-2-aminocyclopentane carboxylic acid", "", "C6H11NO2")
        "1(N)-Methyl-Kynurenine; 1-Methyl-(S)-2-Amino-...; C11H14N2O3"
            → ("1(N)-Methyl-Kynurenine", "1-Methyl-(S)-2-Amino-...", "C11H14N2O3")
        "1-Naphthylalanine, alpha-Naphthylalanine"
            → ("1-Naphthylalanine", "alpha-Naphthylalanine", "")
    """
    if not raw:
        return ("", "", "")

    # Primary split on ";"
    parts = [p.strip() for p in raw.split(";") if p.strip()]

    if not parts:
        return (raw.strip(), "", "")

    # The last token is the formula if it matches the pattern
    formula = ""
    if parts and _looks_like_formula(parts[-1]):
        formula = parts.pop()

    # If only one part remains and it contains a comma, try splitting on commas
    # outside parentheses (handles "Full Name, Synonym" and "Full Name, Formula" patterns)
    if len(parts) == 1 and "," in parts[0]:
        comma_parts = _split_outside_parens(parts[0], ",")
        # Check if last comma-part is a formula
        if comma_parts and not formula and _looks_like_formula(comma_parts[-1]):
            formula = comma_parts.pop()
        if comma_parts:
            parts = comma_parts  # re-assign with comma-split tokens

    full_name = parts[0] if parts else ""
    synonyms  = "; ".join(parts[1:]) if len(parts) > 1 else ""

    return (full_name, synonyms, formula)


# ── API ───────────────────────────────────────────────────────────────────────

def fetch_unusual_aa() -> list[dict]:
    """Fetch all unusual amino acids from the DBAASP API."""
    log.info(f"Fetching unusual amino acids from {API_URL} ...")
    try:
        resp = requests.get(API_URL, params=API_PARAMS, headers=API_HEADERS, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        log.info(f"  Received {len(data)} records.")
        return data
    except requests.exceptions.RequestException as e:
        log.error(f"API request failed: {e}")
        sys.exit(1)
    except ValueError:
        log.error("Failed to parse JSON response.")
        sys.exit(1)


# ── Database ──────────────────────────────────────────────────────────────────

def setup_table(conn: sqlite3.Connection):
    conn.execute(CREATE_TABLE_SQL)
    conn.commit()
    log.info("Table `unusual_amino_acids` ready.")


def upsert_records(conn: sqlite3.Connection, records: list[dict]) -> tuple[int, int]:
    """
    Upsert all records. Returns (inserted, skipped_empty).
    """
    cur = conn.cursor()
    inserted = 0
    skipped  = 0

    for rec in records:
        name_raw = (rec.get("name") or "").strip()
        desc_raw = (rec.get("description") or "").strip()

        if not name_raw:
            skipped += 1
            log.debug(f"Skipping record with empty name: {rec}")
            continue

        full_name, synonyms, formula = parse_description(desc_raw)

        cur.execute(UPSERT_SQL, (name_raw, full_name, synonyms, formula, desc_raw))
        inserted += 1

    conn.commit()
    return inserted, skipped


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Scrape DBAASP unusual amino acids into the SQLite database."
    )
    parser.add_argument(
        "--db", default=DEFAULT_DB,
        help=f"Path to SQLite database (default: {DEFAULT_DB})"
    )
    args = parser.parse_args()

    records = fetch_unusual_aa()

    import os
    os.makedirs(os.path.dirname(args.db) if os.path.dirname(args.db) else ".", exist_ok=True)

    conn = sqlite3.connect(args.db)
    try:
        setup_table(conn)
        inserted, skipped = upsert_records(conn, records)
        log.info(f"Done. Inserted/updated: {inserted} | Skipped (empty name): {skipped}")

        # Quick verification
        count = conn.execute("SELECT COUNT(*) FROM unusual_amino_acids").fetchone()[0]
        log.info(f"Total rows in unusual_amino_acids table: {count}")

        # Show a few examples
        log.info("Sample rows:")
        for row in conn.execute(
            "SELECT name, full_name, molecular_formula FROM unusual_amino_acids "
            "WHERE full_name != '' LIMIT 5"
        ).fetchall():
            log.info(f"  {row[0]!r:25s} → {row[1]!r:40s} [{row[2]}]")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
