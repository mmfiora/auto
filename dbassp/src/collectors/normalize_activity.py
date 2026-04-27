# normalize_activity.py
# Normalizer:
# - Compute peptide MW from SEQUENCE (+ optional N- and C-terminus mods)
# - Split concentration into lower/upper (handles "a-b", ">x", "<x", "±", single value)
# - Convert each bound to µg/mL (primary unit)
# - Classify activity:  active (<=32), not active (>32), unknown (ambiguous >X<=32),
#   or DROPPED (>X where X<=32 — cannot classify)
# - Split targetSpecies into species and strain
# - Write activity_normalized.csv
# + Join por NEW_SEQ con list_min.txt para curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2

import csv
import logging
import re
from src.core.config import Config
from src.core.exceptions import FileProcessingError, DataValidationError

logger = logging.getLogger("dbaasp_pipeline")

# ─── Threshold ────────────────────────────────────────────────────────────────
MIC_THRESHOLD_UGML = 32.0   # µg/mL — active if MIC <= this value


# ─── MW / sequence helpers ────────────────────────────────────────────────────

def calc_mw(seq: str, nterm: str | None, cterm: str | None) -> float:
    """Calculate peptide molecular weight (Da)."""
    mw = sum(Config.AA_MASS.get(a.upper(), 110.0) for a in (seq or "")) + Config.H2O_MASS
    if nterm and nterm.upper() in Config.NTERM_MASS:
        mw += Config.NTERM_MASS[nterm.upper()]
    if cterm and cterm.upper() in Config.CTERM_MASS:
        mw += Config.CTERM_MASS[cterm.upper()]
    return mw


def get_z_prefix(nterm: str | None) -> str:
    """
    Calculate the number of Z's based on N-terminus.
    Each Z represents a C4 block.
    Examples: C4->Z, C8->ZZ, C12->ZZZ, C16->ZZZZ, C20->ZZZZZ
    """
    if not nterm:
        logger.warning("No N-terminus provided, defaulting to C16 (ZZZZ)")
        return "ZZZZ"
    match = re.search(r'C(\d+)', nterm.upper())
    if match:
        num = int(match.group(1))
        z_count = num // 4
        if z_count > 0:
            return "Z" * z_count
        else:
            logger.warning(f"Invalid C value {num}, must be multiple of 4. Defaulting to C16 (ZZZZ)")
            return "ZZZZ"
    logger.warning(f"Could not parse N-terminus '{nterm}', defaulting to C16 (ZZZZ)")
    return "ZZZZ"


# ─── Concentration parsing ─────────────────────────────────────────────────────

def parse_conc(val: str | None, row_num: int | None = None) -> tuple[str, str, int]:
    """
    Split concentration string into (lower_raw, upper_raw, conc_gt).

    conc_gt = 1 means the original value was ">X" or ">=X".

    Semantics:
      - Single  X      → (X,   "",  0)   lower = X, no upper
      - Range   A-B    → (A,   B,   0)
      - Arrow   A->B   → (A,   B,   0)
      - >X / >=X       → (X,   "",  1)   lower = X, flagged as "greater than"
      - <X / <=X       → ("",  X,   0)   only upper bound known
      - mean±err       → (mean-err, mean+err, 0)
    """
    if not val:
        return ("", "", 0)

    s = re.sub(r'\s+', ' ', str(val).strip())

    # Skip malformed values like "4.5.5"
    if s.count('.') > 2 or (s.count('.') > 1 and '±' not in s and '->' not in s and '-' not in s):
        row_info = f" (CSV row {row_num})" if row_num else ""
        logger.warning(f"Malformed concentration value skipped: '{s}'{row_info}")
        return ("", "", 0)

    # mean ± error
    if "±" in s:
        parts = s.split("±", 1)
        try:
            mu  = float(parts[0].strip())
            err = float(parts[1].strip())
            lower = max(0.0, mu - err)
            upper = mu + err
            return (f"{lower:.6g}", f"{upper:.6g}", 0)
        except ValueError:
            row_info = f" (CSV row {row_num})" if row_num else ""
            logger.warning(f"Invalid mean±error format: '{s}'{row_info}")
            return ("", "", 0)

    # Arrow range "A->B"
    if "->" in s:
        a, b = s.split("->", 1)
        return (a.strip(), b.strip(), 0)

    # >=X or >X  →  lower = X, conc_gt = 1
    if s.startswith(">="):
        return (s[2:].strip(), "", 1)
    if s.startswith(">"):
        return (s[1:].strip(), "", 1)

    # <=X or <X  →  upper = X only
    if s.startswith("<="):
        return ("", s[2:].strip(), 0)
    if s.startswith("<"):
        return ("", s[1:].strip(), 0)

    # A-B interval (avoid splitting negative numbers)
    if "-" in s and not s.startswith("-"):
        parts = s.split("-", 1)
        if len(parts) == 2 and parts[0].strip() and parts[1].strip():
            return (parts[0].strip(), parts[1].strip(), 0)

    # Single value → lower = upper = X
    return (s, s, 0)


def to_ugml(val: str, unit: str | None, mw: float | None, row_num: int | None = None) -> str:
    """Convert one numeric concentration value to µg/mL. Returns '' if not possible."""
    if not val or not unit or not mw:
        return ""
    try:
        num = float(val)
    except ValueError:
        row_info = f" (CSV row {row_num})" if row_num else ""
        logger.warning(f"Invalid concentration value: '{val}'{row_info}")
        return ""

    u = (unit or "").strip().lower()
    # Normalise encoding artefacts (µ can arrive as various byte sequences)
    u = u.replace("\xc2\xb5", "µ").replace("\xb5", "µ").replace("?", "µ").replace("\ufffd", "µ")

    # Already µg/mL family — return as-is
    if u in ("µg/ml", "ug/ml", "μg/ml", "mg/l"):
        return f"{num:.6g}"
    if u == "mg/ml":
        return f"{num * 1_000.0:.6g}"
    if u in ("g/l", "g/ml"):
        return f"{num * 1_000.0:.6g}"
    if u in ("ng/ml", "ng/l"):
        return f"{num / 1_000.0:.6g}"

    # Molar units — need MW (Da = g/mol)
    #   µg/mL = µM × (MW g/mol) × (1 µg/µmol) = µM × MW/1000
    if u in ("µm", "um", "μm"):
        return f"{num * (mw / 1_000.0):.6g}"
    if u == "mm":
        return f"{num * 1_000.0 * (mw / 1_000.0):.6g}"
    if u == "nm":
        return f"{num / 1_000.0 * (mw / 1_000.0):.6g}"
    if u == "m":
        return f"{num * 1_000_000.0 * (mw / 1_000.0):.6g}"

    row_info = f" (CSV row {row_num})" if row_num else ""
    logger.warning(f"Unknown unit '{unit}' — cannot convert to µg/mL{row_info}")
    return ""


def classify_activity(lower_ugml: str, upper_ugml: str, conc_gt: int) -> str:
    """
    Classify activity using the 32 µg/mL threshold.

    Rules:
      conc_gt = 0 (exact / range / <X):
        • Use lower_ugml (= the value itself for single; lower end of range).
        • If lower_ugml <= 32  → 'active'
        • If lower_ugml  > 32  → 'not active'
        • If no lower but upper available → use upper as worst case
        • Otherwise → 'unknown'

      conc_gt = 1 (>X):
        • X is the *lower* bound (true MIC is above X).
        • If X  > 32  → 'not active'  (definitively above threshold)
        • If X <= 32  → 'drop'        (cannot classify — MIC could be anywhere > X)
    """
    if conc_gt:
        try:
            x = float(lower_ugml) if lower_ugml else None
        except ValueError:
            x = None
        if x is None:
            return "unknown"
        return "not active" if x > MIC_THRESHOLD_UGML else "drop"

    # Normal case: use lower bound; fall back to upper
    ref = lower_ugml if lower_ugml else upper_ugml
    if not ref:
        return "unknown"
    try:
        v = float(ref)
    except ValueError:
        return "unknown"
    return "active" if v <= MIC_THRESHOLD_UGML else "not active"


# ─── Species splitter ──────────────────────────────────────────────────────────

def split_species(val: str | None) -> tuple[str, str]:
    """
    Split 'Genus species STRAIN' into (species, strain).
    The strain starts at the first token (after word 1) that begins with
    an uppercase letter or a digit.
    """
    if not val:
        return ("", "")
    parts = val.split()
    for i, tok in enumerate(parts[1:], start=1):
        t = tok.strip()
        if t and (t[0].isupper() or t[0].isdigit()):
            return (" ".join(parts[:i]), " ".join(parts[i:]))
    return (val, "")


# ─── Min-list loader ───────────────────────────────────────────────────────────

def load_min_map(path: str = None):
    """
    Read list_min.txt (space/tab-separated) and return:
        {sequence → (curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2)}
    Returns {} on any problem.
    """
    if path is None:
        path = Config.MIN_LIST_FILE

    mapping = {}
    try:
        with open(path, encoding=Config.CSV_ENCODING) as f:
            header = f.readline().strip().split()
            try:
                idx_seq     = header.index("sequence")
                idx_npol    = header.index("npol_min")
                idx_curv    = header.index("curv_min")
                idx_ph      = header.index("pH")
                idx_npol_c0 = header.index("npol_c0")
                idx_npol_c1 = header.index("npol_c1")
                idx_npol_c2 = header.index("npol_c2")
            except ValueError as e:
                logger.warning(f"Missing required columns in {path}: {e}")
                return mapping

            for line in f:
                if not line.strip():
                    continue
                cols = line.strip().split()
                def _get(idx):
                    v = cols[idx] if len(cols) > idx else ""
                    return "" if v == "NA" else v

                seq = _get(idx_seq)
                if seq:
                    mapping[seq] = (
                        _get(idx_curv),
                        _get(idx_npol),
                        _get(idx_ph),
                        _get(idx_npol_c0),
                        _get(idx_npol_c1),
                        _get(idx_npol_c2),
                    )
    except FileNotFoundError:
        logger.warning(f"Min list file not found: {path}")
    except IOError as e:
        logger.error(f"Error reading min list file {path}: {e}")
    return mapping


# ─── Main entry point ──────────────────────────────────────────────────────────

def run(infile: str | None = None, outfile: str | None = None) -> None:
    """
    Read activity.csv, compute MW, split/normalize concentrations to µg/mL,
    classify activity, split species/strain, and write output CSV.

    Rows classified as 'drop' (>X where X <= 32 µg/mL) are excluded from
    the output because the true MIC is unknown relative to the threshold.
    """
    if infile is None:
        infile = Config.OUTPUT_ACTIVITY_CSV
    if outfile is None:
        outfile = Config.OUTPUT_NORMALIZED_CSV

    logger.info(f"Starting activity normalization: {infile} -> {outfile}")

    try:
        min_map = load_min_map()
        logger.info(f"Loaded {len(min_map)} entries from min list")

        with open(infile, encoding=Config.CSV_ENCODING) as f:
            r = csv.DictReader(f)
            fieldnames = list(r.fieldnames) + [
                "MW_Da", "NEW_SEQ",
                "lower_concentration",   # numeric lower bound (raw scale)
                "upper_concentration",   # numeric upper bound (raw scale)
                "conc_gt",               # 1 if original was >X / >=X
                "lower_ugml",            # lower bound in µg/mL
                "upper_ugml",            # upper bound in µg/mL
                "species", "strain",
                "curv_min", "npol_min", "ph_run", "npol_c0", "npol_c1", "npol_c2",
                "activity",              # 'active' | 'not active' | 'unknown'
            ]
            rows = []
            dropped = 0
            row_num = 1

            for row in r:
                row_num += 1
                seq  = (row.get("SEQUENCE") or "").upper()
                n    = row.get("N TERMINUS", "")
                c    = row.get("C TERMINUS", "")
                mw   = calc_mw(seq, n, c)
                unit = row.get("unit", "")

                lo, up, gt = parse_conc(row.get("concentration", ""), row_num)

                row["MW_Da"] = f"{mw:.2f}"
                z_prefix = get_z_prefix(n)
                new_seq  = f"{z_prefix}{seq}{'01' if (c or '').upper() == 'AMD' else '00'}"
                row["NEW_SEQ"] = new_seq

                row["lower_concentration"] = lo
                row["upper_concentration"] = up
                row["conc_gt"]             = gt

                lo_ugml = to_ugml(lo, unit, mw, row_num)
                up_ugml = to_ugml(up, unit, mw, row_num)
                row["lower_ugml"] = lo_ugml
                row["upper_ugml"] = up_ugml

                sp, st = split_species(row.get("targetSpecies", ""))
                row["species"] = sp
                row["strain"]  = st

                # Min-list join — exact NEW_SEQ match only.
                # 00 (AMD) and 01 (non-AMD) are different simulations; never cross-assign.
                curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2 = ("", "", "", "", "", "")
                if new_seq in min_map:
                    curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2 = min_map[new_seq]
                row["curv_min"] = curv_min
                row["npol_min"] = npol_min
                row["ph_run"]   = ph_run
                row["npol_c0"]  = npol_c0
                row["npol_c1"]  = npol_c1
                row["npol_c2"]  = npol_c2

                # Activity classification
                act = classify_activity(lo_ugml, up_ugml, gt)
                if act == "drop":
                    dropped += 1
                    continue   # exclude ambiguous >X rows where X <= 32
                row["activity"] = act
                rows.append(row)

        logger.info(f"Processed {row_num - 1} rows for normalization")
        logger.info(f"Dropped {dropped} rows (>X where X <= {MIC_THRESHOLD_UGML} µg/mL — ambiguous)")
        logger.info(f"Writing {len(rows)} rows to output")

        with open(outfile, "w", encoding=Config.CSV_ENCODING, newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)

        logger.info(f"Successfully wrote normalized data to {outfile}")

    except FileNotFoundError:
        raise FileProcessingError(f"Input file not found: {infile}", filename=infile)
    except csv.Error as e:
        raise FileProcessingError(f"CSV processing error: {e}", filename=infile)
    except IOError as e:
        raise FileProcessingError(f"File I/O error: {e}", filename=outfile)
    except Exception as e:
        logger.error(f"Unexpected error in normalization: {e}")
        raise


if __name__ == "__main__":
    run()
