# physchem.py
import csv
import requests

API_URL = "https://dbaasp.org/peptides/{id}"
HEADERS = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}

def load_ids():
    """
    Read peptide IDs from peptides.csv.
    Accept a column named 'Peptide ID' or 'ID' (case-insensitive).
    Accept values like '51' or 'DBAASPS_51'.
    """
    with open("peptides.csv", encoding="utf-8-sig") as f:
        r = csv.DictReader(f)
        col = None
        for h in r.fieldnames or []:
            if h.lower() in ("peptide id", "id"):
                col = h
                break
        if not col:
            return []
        ids = []
        for row in r:
            raw = (row.get(col) or "").strip()
            if not raw:
                continue
            ids.append(int(raw.split("_")[-1]))
        return ids

def fetch(pid: int):
    """Fetch a single peptide JSON."""
    resp = requests.get(API_URL.format(id=pid), headers=HEADERS, timeout=20)
    resp.raise_for_status()
    return resp.json()

def run():
    ids = load_ids()
    if not ids:
        print("No IDs found in peptides.csv")
        return

    # Fetch all peptide JSONs with progress prints
    data = []
    for i, pid in enumerate(ids, start=1):
        print(f"[{i}/{len(ids)}] Fetching physchem for peptide {pid} ...", end="", flush=True)
        try:
            d = fetch(pid)
            data.append(d)
            print(" ok")
        except Exception as e:
            print(f" fail ({e})")

    if not data:
        print("No data fetched; nothing to write.")
        return

    # Collect property names preserving first-seen order, skip any property titled exactly "ID"
    props, seen = [], set()
    for d in data:
        for p in d.get("physicoChemicalProperties") or []:
            name = (p.get("name") or "").strip()
            if not name or name.upper() == "ID":
                continue
            if name not in seen:
                seen.add(name)
                props.append(name)

    # Final header: base peptide columns + dynamic physchem properties
    header = ["Peptide ID", "N TERMINUS", "SEQUENCE", "C TERMINUS"] + props

    # Write CSV, always one row per peptide (fill missing properties with empty string)
    with open("physchem.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()
        for d in data:
            row = {
                "Peptide ID": str(d.get("id", "")),                          # numeric peptide ID
                "N TERMINUS": (d.get("nTerminus") or {}).get("name", ""),    # only .name
                "SEQUENCE": d.get("sequence", ""),                            # original casing
                "C TERMINUS": (d.get("cTerminus") or {}).get("name", ""),    # only .name
            }
            # Map this peptide's properties
            values = {}
            for p in d.get("physicoChemicalProperties") or []:
                name = (p.get("name") or "").strip()
                if not name or name.upper() == "ID":
                    continue
                values[name] = str(p.get("value", "")).strip()
            # Fill all columns consistently
            for name in props:
                row[name] = values.get(name, "")
            w.writerow(row)


