"""
migrate_db.py
-------------
One-shot migration script for dbassp/data/output/dbaasp_data.sqlite.

Changes applied (v1):
  1. Rename `activities`       → `normalized_activity`
       - Add c_min_uM, c_max_uM (normalized concentration bounds in µM)
       - Add species, strain    (split from target_species)
       - Keep original columns intact
  2. Rename `physchem_properties` → `physchem_properties_aa`
  3. Merge `lipophilicity` data (smiles, logp, logd) into `peptides` table,
     then drop `lipophilicity`
  4. Add molecular_weight (Da) and total_charge columns to `peptides`

Changes applied (v2):
  5. Add to normalized_activity:
       - c_min_ugml  REAL   — min concentration in µg/ml
       - c_max_ugml  REAL   — max concentration in µg/ml (NULL for single values)
       - conc_gt     INTEGER — 1 if original string started with '>' / '>='
       - activity    TEXT   — 'active' / 'not active' / 'unknown'

Concentration column semantics (both µM and µg/ml columns):
  - Single number  X      → c_min = X,       c_max = NULL,      conc_gt = 0
  - Interval A-B / A->B  → c_min = A,        c_max = B,         conc_gt = 0
  - >X or >=X            → c_min = X,        c_max = NULL,      conc_gt = 1
  - <X or <=X            → c_min = NULL,     c_max = X,         conc_gt = 0
  - ± mean/error         → c_min = mean-err, c_max = mean+err,  conc_gt = 0

Activity threshold: 32 µg/ml (applied to c_min_ugml)
  - conc_gt=0: active if <= 32, not active if > 32
  - conc_gt=1: unknown if <= 32, not active if > 32
  - c_min_ugml NULL: unknown

Units supported: µM/uM, mM, nM, M, µg/ml/ug/ml/mg/L, g/L, ng/mL
"""

import re
import sqlite3
import logging

DB_PATH = "dbassp/data/output/dbaasp.sqlite"
logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
log = logging.getLogger("migrate_db")

# ─── Amino-acid masses (monoisotopic-like, as used in the pipeline) ────────────
AA_MASS = {
    "A": 71.08,  "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.15,
    "E": 129.12, "Q": 128.13, "G": 57.05,  "H": 137.14, "I": 113.16,
    "L": 113.16, "K": 128.17, "M": 131.20, "F": 147.18, "P": 97.12,
    "S": 87.08,  "T": 101.11, "W": 186.21, "Y": 163.18, "V": 99.13,
    "Z": 56.10,  "X": 110.0,
}
H2O_MASS   = 18.02
NTERM_MASS = {"C16": 239.2, "C14": 211.2, "C12": 183.2, "C10": 155.2,
               "C8": 127.2, "C18": 267.2, "C20": 295.2, "C6": 99.2,
               "C4": 71.2}
CTERM_MASS = {"AMD": -0.98}

# Amino-acid net charge at pH 7 (simplified: K,R,H=+1 ; D,E=-1)
AA_CHARGE = {
    "K": 1, "R": 1, "H": 1,
    "D": -1, "E": -1,
}
NTERM_CHARGE = {}   # termini don't contribute extra charge in this model
CTERM_CHARGE = {}


def calc_mw(seq: str, nterm: str, cterm: str) -> float:
    mw = sum(AA_MASS.get(a.upper(), 110.0) for a in (seq or "")) + H2O_MASS
    key_n = (nterm or "").upper()
    key_c = (cterm or "").upper()
    mw += NTERM_MASS.get(key_n, 0.0)
    mw += CTERM_MASS.get(key_c, 0.0)
    return round(mw, 4)


def calc_charge(seq: str) -> int:
    return sum(AA_CHARGE.get(a.upper(), 0) for a in (seq or ""))


# ─── Concentration parsing ─────────────────────────────────────────────────────

def parse_conc(val: str) -> tuple[str, str]:
    """Return (lower_raw, upper_raw) numeric strings (empty string = not set)."""
    if not val:
        return ("", "")
    s = re.sub(r"\s+", " ", str(val).strip())

    # mean ± error  →  (mean-err, mean+err)
    if "±" in s:
        parts = s.split("±", 1)
        try:
            mu  = float(parts[0].strip())
            err = float(parts[1].strip())
            return (f"{max(0.0, mu-err):.6g}", f"{mu+err:.6g}")
        except ValueError:
            return ("", "")

    # arrow range "A->B"
    if "->" in s:
        a, b = s.split("->", 1)
        return (a.strip(), b.strip())

    # >=X  →  (None, X)
    if s.startswith(">="):
        return ("", s[2:].strip())

    # >X  →  (None, X)   i.e. true MIC is above X  →  c_max = X (as lower bound of "above")
    if s.startswith(">"):
        return ("", s[1:].strip())

    # <=X or <X  →  (None, X)
    if s.startswith("<="):
        return ("", s[2:].strip())
    if s.startswith("<"):
        return ("", s[1:].strip())

    # A-B range (skip negative numbers / leading dash)
    if "-" in s and not s.startswith("-"):
        parts = s.split("-", 1)
        if len(parts) == 2 and parts[0].strip() and parts[1].strip():
            return (parts[0].strip(), parts[1].strip())

    # Single value
    return (s, s)


def to_uM(val: str, unit: str, mw: float | None) -> float | None:
    """Convert a numeric concentration string to µM. Returns None on failure."""
    if not val:
        return None
    try:
        num = float(val)
    except ValueError:
        return None

    u = (unit or "").strip().lower()
    # Normalise encoding artefacts:
    #  µ can be stored as \xb5, \xc2\xb5, or the Unicode replacement char \ufffd
    u = u.replace("\xc2\xb5", "µ").replace("\xb5", "µ").replace("?", "µ")
    u = u.replace("\ufffd", "µ")   # replacement char from broken UTF-8 round-trip

    if u in ("µm", "um", "μm"):
        return round(num, 6)
    if u == "mm":
        return round(num * 1_000.0, 6)
    if u == "nm":
        return round(num / 1_000.0, 6)
    if u == "m":
        return round(num * 1_000_000.0, 6)
    if u in ("µg/ml", "ug/ml", "μg/ml", "mg/l", "mg/ml"):
        if not mw:
            return None
        return round(num / (mw / 1_000.0), 6)
    if u in ("g/l", "g/ml"):
        if not mw:
            return None
        return round((num * 1_000.0) / (mw / 1_000.0), 6)
    if u in ("ng/ml", "ng/l"):
        if not mw:
            return None
        return round((num / 1_000.0) / (mw / 1_000.0), 6)
    # unknown unit
    log.warning(f"Unknown unit '{unit}' — cannot convert to µM")
    return None


def split_species(val: str) -> tuple[str, str]:
    """
    Split 'Genus species STRAIN' into (species, strain).
    Rule: the strain starts at the first token (after the 1st word) that
    begins with an uppercase letter or a digit.
    Examples:
      'Escherichia coli ATCC 25922'   → ('Escherichia coli', 'ATCC 25922')
      'Candida albicans TIMM 0144'    → ('Candida albicans', 'TIMM 0144')
      'Candida krusei'                → ('Candida krusei', '')
    """
    if not val:
        return ("", "")
    parts = val.split()
    for i, tok in enumerate(parts[1:], start=1):
        t = tok.strip()
        if t and (t[0].isupper() or t[0].isdigit()):
            return (" ".join(parts[:i]), " ".join(parts[i:]))
    return (val, "")


# ─── Migration ─────────────────────────────────────────────────────────────────

def run_migration():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = OFF")
    cur  = conn.cursor()

    # ── Add molecular_weight and total_charge to peptides ─────────────────────
    log.info("Adding molecular_weight / total_charge to peptides …")
    existing_cols = {r[1] for r in cur.execute("PRAGMA table_info(peptides)").fetchall()}

    for col, typ in [("molecular_weight", "REAL"), ("total_charge", "INTEGER")]:
        if col not in existing_cols:
            cur.execute(f"ALTER TABLE peptides ADD COLUMN {col} {typ}")
            log.info(f"  Added column: {col}")

    cur.execute("SELECT id, sequence FROM peptides")
    rows = cur.fetchall()
    updated = 0
    for pid, seq in rows:
        mw  = calc_mw(seq or "", "", "")   # termini already factored into peptides.n_terminus / c_terminus
        # Re-read termini for accurate MW
        cur2 = conn.cursor()
        cur2.execute("SELECT n_terminus, c_terminus FROM peptides WHERE id=?", (pid,))
        r = cur2.fetchone()
        nterm, cterm = (r[0] or "", r[1] or "") if r else ("", "")
        mw  = calc_mw(seq or "", nterm, cterm)
        chg = calc_charge(seq or "")
        cur.execute(
            "UPDATE peptides SET molecular_weight=?, total_charge=? WHERE id=?",
            (mw, chg, pid)
        )
        updated += 1

    conn.commit()
    conn.execute("PRAGMA foreign_keys = ON")
    conn.close()
    log.info(f"  Updated {updated} peptides with MW and charge")
    log.info("Migration complete ✓")



# ─── µg/ml conversion helper ───────────────────────────────────────────────────

def _to_ugml(val: str, unit: str, mw: float | None) -> float | None:
    """Convert a numeric concentration string to µg/ml. Returns None on failure."""
    if not val:
        return None
    try:
        num = float(val)
    except ValueError:
        return None

    u = (unit or "").strip().lower()
    u = u.replace("\xc2\xb5", "µ").replace("\xb5", "µ").replace("?", "µ")
    u = u.replace("\ufffd", "µ")

    if u in ("µg/ml", "ug/ml", "μg/ml", "mg/l"):
        return round(num, 6)
    if u == "mg/ml":
        return round(num * 1_000.0, 6)
    if u in ("g/l", "g/ml"):
        return round(num * 1_000.0, 6)
    if u in ("ng/ml", "ng/l"):
        return round(num / 1_000.0, 6)
    # molar units — need MW to convert
    if not mw:
        return None
    if u in ("µm", "um", "μm"):
        return round(num * (mw / 1_000.0), 6)
    if u == "mm":
        return round(num * 1_000.0 * (mw / 1_000.0), 6)
    if u == "nm":
        return round(num / 1_000.0 * (mw / 1_000.0), 6)
    if u == "m":
        return round(num * 1_000_000.0 * (mw / 1_000.0), 6)
    return None


def _classify_activity(c_min_ugml: float | None, conc_gt: int) -> str:
    """Return 'active', 'not active', or 'unknown' based on µg/ml threshold of 32."""
    if c_min_ugml is None:
        return "unknown"
    if conc_gt:
        # form was ">X": if X > 32 → definitely not active; if X <= 32 → unknown
        return "not active" if c_min_ugml > 32 else "unknown"
    else:
        return "active" if c_min_ugml <= 32 else "not active"


def _parse_conc_v2(val: str) -> tuple[str, str, int]:
    """
    Extended parse_conc that also returns conc_gt flag.
    Returns (lower_raw, upper_raw, conc_gt)
      lower_raw / upper_raw: numeric strings (empty = not set)
      conc_gt: 1 if form was '>X' or '>=X', else 0
    """
    if not val:
        return ("", "", 0)
    s = re.sub(r"\s+", " ", str(val).strip())

    # mean ± error
    if "±" in s:
        parts = s.split("±", 1)
        try:
            mu  = float(parts[0].strip())
            err = float(parts[1].strip())
            return (f"{max(0.0, mu-err):.6g}", f"{mu+err:.6g}", 0)
        except ValueError:
            return ("", "", 0)

    # arrow range
    if "->" in s:
        a, b = s.split("->", 1)
        return (a.strip(), b.strip(), 0)

    # >=X or >X  → c_min = X, c_max = NULL, conc_gt = 1
    if s.startswith(">="):
        return (s[2:].strip(), "", 1)
    if s.startswith(">"):
        return (s[1:].strip(), "", 1)

    # <=X or <X  → c_min = NULL, c_max = X
    if s.startswith("<="):
        return ("", s[2:].strip(), 0)
    if s.startswith("<"):
        return ("", s[1:].strip(), 0)

    # A-B interval
    if "-" in s and not s.startswith("-"):
        parts = s.split("-", 1)
        if len(parts) == 2 and parts[0].strip() and parts[1].strip():
            return (parts[0].strip(), parts[1].strip(), 0)

    # Single value → c_min = X, c_max = NULL
    return (s, "", 0)


# ─── Migration v2 ──────────────────────────────────────────────────────────────

def migrate_v2():
    """
    Adds c_min_ugml, c_max_ugml, conc_gt, activity to normalized_activity
    and back-fills them for every existing row.
    """
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = OFF")
    cur  = conn.cursor()

    # ── Add new columns (idempotent) ───────────────────────────────────────────
    existing_cols = {r[1] for r in cur.execute("PRAGMA table_info(normalized_activity)").fetchall()}
    new_cols = [
        ("c_min_ugml", "REAL"),
        ("c_max_ugml", "REAL"),
        ("conc_gt",    "INTEGER"),
        ("activity",   "TEXT"),
    ]
    for col, typ in new_cols:
        if col not in existing_cols:
            cur.execute(f"ALTER TABLE normalized_activity ADD COLUMN {col} {typ}")
            log.info(f"  Added column: {col}")

    # ── Back-fill ──────────────────────────────────────────────────────────────
    cur.execute("""
        SELECT na.id, na.concentration, na.unit, p.molecular_weight
        FROM normalized_activity na
        JOIN peptides p ON p.id = na.peptide_id
    """)
    rows = cur.fetchall()
    updated = 0
    for row_id, conc, unit, mw in rows:
        lo, hi, gt = _parse_conc_v2(conc or "")
        c_min_ugml = _to_ugml(lo, unit, mw)
        c_max_ugml = _to_ugml(hi, unit, mw)
        activity   = _classify_activity(c_min_ugml, gt)
        cur.execute("""
            UPDATE normalized_activity
            SET c_min_ugml=?, c_max_ugml=?, conc_gt=?, activity=?
            WHERE id=?
        """, (c_min_ugml, c_max_ugml, gt, activity, row_id))
        updated += 1

    conn.commit()
    conn.execute("PRAGMA foreign_keys = ON")
    conn.close()
    log.info(f"  Updated {updated} normalized_activity rows with µg/ml + activity")
    log.info("Migration v2 complete ✓")


if __name__ == "__main__":
    run_migration()
    migrate_v2()
