# normalize_activity.py
# Minimal normalizer:
# - Compute peptide MW from SEQUENCE (+ optional N- and C-terminus mods)
# - Split concentration into lower/upper (handles "a-b", ">x", "<x", single value)
# - Convert each bound to µM
# - Split targetSpecies into species and strain using a simple heuristic
# - Write activity_normalized.csv (or a custom outfile)
# + Join por NEW_SEQ con list_min.txt para curv_min y npol_min

import csv
import re
from config import Config

def calc_mw(seq: str, nterm: str | None, cterm: str | None) -> float:
    """Calculate peptide molecular weight (Da)."""
    mw = sum(Config.AA_MASS.get(a.upper(), 110.0) for a in (seq or "")) + Config.H2O_MASS
    if nterm and nterm.upper() in Config.NTERM_MASS:
        mw += Config.NTERM_MASS[nterm.upper()]
    if cterm and cterm.upper() in Config.CTERM_MASS:
        mw += Config.CTERM_MASS[cterm.upper()]
    return mw

def parse_conc(val: str | None) -> tuple[str, str]:
    """Split concentration string into (lower, upper) strings."""
    if not val:
        return ("", "")
    s = val.strip()
    if "-" in s:
        a, b = s.split("-", 1)
        return (a.strip(), b.strip())
    if s.startswith(">"):
        return (s[1:].strip(), "")
    if s.startswith("<"):
        return ("", s[1:].strip())
    return (s, s)

def to_uM(val: str, unit: str | None, mw: float | None) -> str:
    """Convert one numeric concentration value to µM. Returns '' if not possible."""
    if not val or not unit or not mw:
        return ""
    try:
        num = float(val)
    except Exception:
        return ""
    u = unit.lower()
    if u in ("µm", "um"):
        return f"{num:.6g}"
    if u == "mm":
        return f"{(num * 1000.0):.6g}"
    if u == "nm":
        return f"{(num / 1000.0):.6g}"
    if u == "m":
        return f"{(num * 1_000_000.0):.6g}"
    # mass/volume units require MW (Da ~ g/mol)
    if u in ("µg/ml", "ug/ml", "mg/l"):
        return f"{(num / (mw / 1000.0)):.6g}"
    if u == "g/l":
        return f"{((num * 1000.0) / (mw / 1000.0)):.6g}"
    if u == "ng/ml":
        return f"{((num / 1000.0) / (mw / 1000.0)):.6g}"
    return ""

def split_species(val: str | None) -> tuple[str, str]:
    """
    Split target species into species and strain with a simple rule:
    - skip the first word (usually the genus)
    - if any later token starts with an uppercase letter OR a digit, that token begins the strain
    """
    if not val:
        return ("", "")
    parts = val.split()
    for i, tok in enumerate(parts[1:], start=1):  # skip first word
        t = tok.strip()
        if not t:
            continue
        if t[0].isupper() or t[0].isdigit():
            return (" ".join(parts[:i]), " ".join(parts[i:]))
    return (val, "")

# --- NUEVO: cargar mapa {NEW_SEQ -> (curv_min, npol_min)} desde list_min.txt ---
def load_min_map(path: str = None):
    """
    Lee list_min.txt (separado por espacios/tabs) y devuelve un dict:
        key = sequence (por ej., 'ZZZZKLK01')
        value = (curv_min, npol_min) como strings ('' si NA)
    Si hay problemas, devuelve {}.
    """
    if path is None:
        path = Config.MIN_LIST_FILE
        
    mapping = {}
    try:
        with open(path, encoding=Config.CSV_ENCODING) as f:
            header = f.readline().strip().split()
            idx_seq = header.index("sequence")
            idx_npol = header.index("npol_min")
            idx_curv = header.index("curv_min")
            for line in f:
                if not line.strip():
                    continue
                cols = line.strip().split()
                # robustez ante 'NA'
                seq = cols[idx_seq] if len(cols) > idx_seq else ""
                npol = cols[idx_npol] if len(cols) > idx_npol else ""
                curv = cols[idx_curv] if len(cols) > idx_curv else ""
                if npol == "NA":
                    npol = ""
                if curv == "NA":
                    curv = ""
                if seq:
                    mapping[seq] = (curv, npol)
    except Exception:
        pass
    return mapping

def run(infile: str = None, outfile: str = None) -> None:
    """Read activity.csv, compute MW, split/normalize concentrations, split species/strain, write output CSV."""
    if infile is None:
        infile = Config.OUTPUT_ACTIVITY_CSV
    if outfile is None:
        outfile = Config.OUTPUT_NORMALIZED_CSV
        
    # NUEVO: cargar join por NEW_SEQ
    min_map = load_min_map()

    with open(infile, encoding=Config.CSV_ENCODING) as f:
        r = csv.DictReader(f)
        fieldnames = list(r.fieldnames) + [
            "MW_Da", "NEW_SEQ",
            "lower_concentration", "upper_concentration",
            "lower_uM", "upper_uM",
            "species", "strain",
            # --- NUEVAS columnas ---
            "curv_min", "npol_min",
        ]
        rows = []
        for row in r:
            seq = (row.get("SEQUENCE") or "").upper()
            n = row.get("N TERMINUS", "")
            c = row.get("C TERMINUS", "")
            mw = calc_mw(seq, n, c)

            lo, up = parse_conc(row.get("concentration", ""))

            row["MW_Da"] = f"{mw:.2f}"
            new_seq = f"ZZZZ{seq}{'01' if (c or '').upper() == 'AMD' else '00'}"
            row["NEW_SEQ"] = new_seq
            row["lower_concentration"] = lo
            row["upper_concentration"] = up
            row["lower_uM"] = to_uM(lo, row.get("unit", ""), mw)
            row["upper_uM"] = to_uM(up, row.get("unit", ""), mw)

            sp, st = split_species(row.get("targetSpecies", ""))
            row["species"] = sp
            row["strain"] = st

            # --- Join con curv_min / npol_min por NEW_SEQ ---
            curv_min, npol_min = ("", "")
            if new_seq in min_map:
                curv_min, npol_min = min_map[new_seq]
            row["curv_min"] = curv_min
            row["npol_min"] = npol_min

            rows.append(row)

    with open(outfile, "w", encoding=Config.CSV_ENCODING, newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    return  # explicit return to avoid accidental fall-through

#if __name__ == "__main__":
    run()
