# extract_physchem.py
import re
import csv
import os
from typing import Dict, List, Optional

PROPS = [
    "Normalized Hydrophobicity",
    "Net Charge",
    "Isoelectric Point",
    "Penetration Depth",
    "Tilt Angle",
    "Disordered Conformation Propensity",
    "Linear Moment",
    "Propensity to in vitro Aggregation",
    "Angle Subtended by the Hydrophobic Residues",
    "Amphiphilicity Index",
    "Propensity to PPII coil",
]

NUM_RE = re.compile(r"[-+]?\d+(?:[.,]\d+)?")

def _to_float(s: str) -> float:
    return float(s.replace(",", "."))

def extract_physchem(filename: str) -> Optional[Dict[str, float]]:
    """
    Simple extractor: take numbers from the Physico-Chemical Properties block
    and map them in order to PROPS.
    """
    with open(filename, encoding="utf-8") as f:
        text = f.read()

    # Start of section
    m_start = re.search(r"Physico-?Chemical Properties", text, re.I)
    if not m_start:
        print("No se encontró la sección Physico-Chemical Properties")
        return None

    start = m_start.end()

    # End of section: next common header or end of file
    m_end = re.search(
        r"\n(?:Activity against|General Information|Amino acid composition|Structure|References|$)",
        text[start:], re.I
    )
    end = start + (m_end.start() if m_end else len(text) - start)

    block = text[start:end]

    # Extract all numbers in order
    nums = NUM_RE.findall(block)
    if len(nums) < len(PROPS):
        print(f"No hay suficientes valores: {len(nums)} encontrados, {len(PROPS)} requeridos")
        return None

    values = list(map(_to_float, nums[:len(PROPS)]))
    return dict(zip(PROPS, values))

def write_physchem_csv(txt_files: List[str], csv_path: str) -> int:
    """
    Write a single CSV with header ['Peptide ID'] + PROPS.
    Returns number of rows written.
    """
    header = ["Peptide ID"] + PROPS
    written = 0
    # Ensure parent directory exists (optional)
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)

    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for path in txt_files:
            data = extract_physchem(path)
            if not data:
                continue
            base = os.path.basename(path)
            m = re.search(r"peptide_(\d+)\.txt$", base, re.I)
            pid = m.group(1) if m else ""
            row = {"Peptide ID": pid}
            row.update(data)
            w.writerow(row)
            written += 1
    return written

# Test rápido opcional
#if __name__ == "__main__":
    files = ["peptide_51.txt"]
    n = write_physchem_csv(files, "physchem.csv")
    print(f"Filas escritas: {n}")

