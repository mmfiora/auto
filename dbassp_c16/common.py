# common.py
# Shared functions and constants for DBAASP API interaction

import csv
import requests

API_URL = "https://dbaasp.org/peptides/{id}"
HEADERS = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}

def load_ids(csv_file: str = "peptides.csv"):
    """
    Read peptide IDs from CSV file.
    Accept a column named 'Peptide ID' or 'ID' (case-insensitive).
    Accept values like '51' or 'DBAASPS_51'.
    """
    with open(csv_file, encoding="utf-8-sig") as f:
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
    """Fetch a single peptide JSON from DBAASP API."""
    resp = requests.get(API_URL.format(id=pid), headers=HEADERS, timeout=20)
    resp.raise_for_status()
    return resp.json()