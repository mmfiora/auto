# extract_activity.py
import os
import re
import csv
from typing import List, Dict, Tuple, Optional, Iterable

HEADERS = [
    "Peptide ID",
    "N TERMINUS",
    "SEQUENCE",
    "C TERMINUS",
    "Target Species",
    "Activity Measure",
    "Activity",
    "Unit",
    "pH",
    "Ionic Strength mM",
    "Salt Type",
    "Medium",
    "CFU",
    "Note",
    "Reference",
]

# ---- Activity parsing (unchanged) ----
MEASURE_RE = re.compile(r"\b(MIC|MBC|IC50|MBEC|LD50|ED50|Hemolysis)\b", re.I)
UNIT_RE = re.compile(r"(?:μ|µ|u)?[pnm]?M|(?:μ|µ|u)?g/?mL|ng/?mL|mg/?L|%|ppm", re.I)
CFU_RE = re.compile(r"\b(?:\d+(?:[Ee]\d+)|10\^\d+)\b")
PH_RE = re.compile(r"\bpH\s*(\d+(?:\.\d+)?)", re.I)
MM_RE = re.compile(r"(\d+(?:\.\d+)?)\s*mM", re.I)
RANGE_NUM_RE = re.compile(r"[<>]?\s*\d+(?:\.\d+)?(?:\s*-\s*\d+(?:\.\d+)?)?")
MEDIUM_RE = re.compile(r"\b([A-Z]{2,}(?:-\d{1,4})?)\b")
NOT_MEDIUM = {"MIC", "MBC", "IC50", "MBEC", "LD50", "ED50", "CFU", "ATCC"}

def _clean(s: str) -> str:
    return " ".join(s.split())

def parse_line(line: str) -> Dict[str, str]:
    """Parse a single activity line into the fixed schema. Drop hemolysis rows."""
    row = {
        "Target Species": "",
        "Activity Measure": "",
        "Activity": "",
        "Unit": "",
        "pH": "",
        "Ionic Strength mM": "",
        "Salt Type": "",
        "Medium": "",
        "CFU": "",
        "Note": "",
        "Reference": "",
    }
    text = line.strip()
    if not text:
        return {}

    m_meas = MEASURE_RE.search(text)
    if not m_meas:
        return {}
    measure = m_meas.group(1).capitalize()
    if measure.lower() == "hemolysis":
        return {}

    row["Activity Measure"] = measure
    row["Target Species"] = text[:m_meas.start()].strip()

    tail = text[m_meas.end():].strip()

    m_num = RANGE_NUM_RE.search(tail)
    if m_num:
        row["Activity"] = m_num.group(0).replace(" ", "")
        tail = tail[m_num.end():].strip()

    m_unit = UNIT_RE.search(tail)
    if m_unit:
        row["Unit"] = m_unit.group(0)
        tail = tail[m_unit.end():].strip()

    m_ph = PH_RE.search(text)
    if m_ph:
        row["pH"] = m_ph.group(1)
    m_mM = MM_RE.search(text)
    if m_mM:
        row["Ionic Strength mM"] = m_mM.group(1)

    m_ref_final = re.search(r"(\d+)\s*$", tail)
    ref = m_ref_final.group(1) if m_ref_final else ""

    medium = ""
    note_work = tail
    m_med = MEDIUM_RE.search(note_work)
    if m_med:
        token = m_med.group(1)
        if token not in NOT_MEDIUM:
            medium = token
            note_work = (note_work[:m_med.start()] + note_work[m_med.end():]).strip()

    m_cfu = CFU_RE.search(note_work)
    if m_cfu:
        row["CFU"] = m_cfu.group(0)
        note_work = (note_work[:m_cfu.start()] + note_work[m_cfu.end():]).strip()

    m_ref = re.search(r"(\d+)\s*$", note_work)
    if m_ref:
        row["Reference"] = m_ref.group(1)
        note_work = note_work[:m_ref.start()].strip()
    else:
        row["Reference"] = ref

    row["Medium"] = medium
    row["Note"] = _clean(note_work)
    return row

def extract_activity_from_text(text: str) -> List[Dict[str, str]]:
    """Extract and parse the 'Activity against target species' block from full page text."""
    m_start = re.search(r"Activity against target species", text, re.I)
    if not m_start:
        return []
    start = m_start.end()
    m_end = re.search(
        r"\n(?:Physico-?Chemical Properties|General Information|Amino acid composition|Structure|References|$)",
        text[start:], re.I
    )
    end = start + (m_end.start() if m_end else len(text) - start)
    block = text[start:end].strip()

    lines = [ln.strip() for ln in block.splitlines() if ln.strip()]
    if not lines:
        return []

    cleaned = [ln for ln in lines if ln not in HEADERS]

    rows: List[Dict[str, str]] = []
    for ln in cleaned:
        r = parse_line(ln)
        if r and r.get("Activity Measure", "").lower() != "hemolysis":
            rows.append(r)
    return rows

def extract_activity(filename: str) -> List[Dict[str, str]]:
    """Return activity rows for a given peptide_<id>.txt file."""
    with open(filename, encoding="utf-8") as f:
        text = f.read()
    return extract_activity_from_text(text)

# ---- Basic info from peptides.csv (source of truth) ----
def _norm(s: Optional[str]) -> str:
    return (s or "").strip()

def _norm_key(s: Optional[str]) -> str:
    return _norm(s).casefold()

def _best_delimiter(sample: str) -> str:
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        return dialect.delimiter
    except Exception:
        return ","  # default

def load_basic_info_from_csv(csv_path: str) -> Dict[str, Dict[str, str]]:
    """
    Load mapping from peptides.csv keyed by Peptide ID (string of digits).
    Accepts delimiters: comma, semicolon, tab. Header names are normalized case-insensitively.
    Expected logical columns: 'Peptide ID', 'N TERMINUS', 'SEQUENCE', 'C TERMINUS'.
    """
    lookup: Dict[str, Dict[str, str]] = {}
    if not csv_path or not os.path.isfile(csv_path):
        return lookup

    with open(csv_path, "r", encoding="utf-8-sig", newline="") as fh:
        head = fh.read(2048)
        fh.seek(0)
        delim = _best_delimiter(head)
        reader = csv.DictReader(fh, delimiter=delim)

        # Build a case-insensitive header map
        header_map = { _norm_key(h): h for h in (reader.fieldnames or []) }

        def col(*names: Iterable[str]) -> Optional[str]:
            for n in names:
                if _norm_key(n) in header_map:
                    return header_map[_norm_key(n)]
            return None

        col_id = col("Peptide ID", "PeptideID", "ID")
        col_n  = col("N TERMINUS", "NTerminus", "N_Terminus", "N-TERMINUS")
        col_s  = col("SEQUENCE", "Sequence")
        col_c  = col("C TERMINUS", "CTerminus", "C_Terminus", "C-TERMINUS")

        if not all([col_id, col_n, col_s, col_c]):
            # Headers missing; return empty so caller leaves fields blank
            return {}

        for row in reader:
            pid = _norm(row.get(col_id))
            if not pid:
                continue
            lookup[pid] = {
                "N TERMINUS": _norm(row.get(col_n)),
                "SEQUENCE": _norm(row.get(col_s)),
                "C TERMINUS": _norm(row.get(col_c)),
            }
    return lookup

# ---- Writer ----
def write_activity_csv(
    txt_files: List[str],
    csv_path: str,
    peptides_csv: Optional[str] = "peptides.csv",
) -> int:
    """
    Write one CSV with one header row and all parsed activity rows from txt_files.
    Adds:
      - 'Peptide ID' from filename 'peptide_<id>.txt'
      - 'N TERMINUS', 'SEQUENCE', 'C TERMINUS' from peptides.csv (no text fallback)
    """
    lookup = load_basic_info_from_csv(peptides_csv) if peptides_csv else {}

    total = 0
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=HEADERS)
        writer.writeheader()
        for txt in txt_files:
            base = os.path.basename(txt)
            m = re.search(r"peptide_(\d+)\.txt$", base, re.I)
            pid = m.group(1) if m else ""

            rows = extract_activity(txt)

            basic = lookup.get(pid, {"N TERMINUS": "", "SEQUENCE": "", "C TERMINUS": ""})

            for r in rows:
                out = {
                    "Peptide ID": pid,
                    "N TERMINUS": basic.get("N TERMINUS", ""),
                    "SEQUENCE": basic.get("SEQUENCE", ""),
                    "C TERMINUS": basic.get("C TERMINUS", ""),
                }
                for h in HEADERS:
                    if h in out:
                        continue
                    out[h] = r.get(h, "")
                writer.writerow(out)
                total += 1
    return total

#if __name__ == "__main__":
    files = ["peptide_51.txt"]
    n = write_activity_csv(files, "activity.csv", peptides_csv="peptides.csv")
    print(f"Wrote {n} rows to activity.csv")

