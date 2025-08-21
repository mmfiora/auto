# normalize_activity.py
# Minimal normalizer:
# - Compute peptide MW from SEQUENCE (+ optional N- and C-terminus mods)
# - Split concentration into lower/upper (handles "a-b", ">x", "<x", single value)
# - Convert each bound to µM
# - Split targetSpecies into species and strain using a simple heuristic
# - Write activity_normalized.csv (or a custom outfile)

import csv
import re

# Average masses (Da)
AA_MASS = {
    "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.15,
    "E": 129.12, "Q": 128.13, "G": 57.05,  "H": 137.14, "I": 113.16,
    "L": 113.16, "K": 128.17, "M": 131.20, "F": 147.18, "P": 97.12,
    "S": 87.08,  "T": 101.11, "W": 186.21, "Y": 163.18, "V": 99.13,
    "Z": 56.10,  # special C4 block
    "X": 110.0   # generic unknown residue
}

# Constants (Da)
H2O = 18.02
NTERM = {"C16": 239.2}     # simple N-terminus additions
CTERM = {"AMD": -0.98}     # small delta vs. free COOH

def calc_mw(seq: str, nterm: str | None, cterm: str | None) -> float:
    """Calculate peptide molecular weight (Da)."""
    mw = sum(AA_MASS.get(a.upper(), 110.0) for a in (seq or "")) + H2O
    if nterm and nterm.upper() in NTERM:
        mw += NTERM[nterm.upper()]
    if cterm and cterm.upper() in CTERM:
        mw += CTERM[cterm.upper()]
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

def run(infile: str = "activity.csv", outfile: str = "activity_normalized.csv") -> None:
    """Read activity.csv, compute MW, split/normalize concentrations, split species/strain, write output CSV."""
    with open(infile, encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        fieldnames = list(r.fieldnames) + [
            "MW_Da", "NEW_SEQ",
            "lower_concentration", "upper_concentration",
            "lower_uM", "upper_uM",
            "species", "strain"
        ]
        rows = []
        for row in r:
            seq = (row.get("SEQUENCE") or "").upper()
            n = row.get("N TERMINUS", "")
            c = row.get("C TERMINUS", "")
            mw = calc_mw(seq, n, c)

            lo, up = parse_conc(row.get("concentration", ""))

            row["MW_Da"] = f"{mw:.2f}"
            row["NEW_SEQ"] = f"ZZZZ{seq}{'01' if (c or '').upper() == 'AMD' else '00'}"
            row["lower_concentration"] = lo
            row["upper_concentration"] = up
            row["lower_uM"] = to_uM(lo, row.get("unit", ""), mw)
            row["upper_uM"] = to_uM(up, row.get("unit", ""), mw)

            sp, st = split_species(row.get("targetSpecies", ""))
            row["species"] = sp
            row["strain"] = st

            rows.append(row)

    with open(outfile, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    return  # explicit return to avoid accidental fall-through

if __name__ == "__main__":
    run()
