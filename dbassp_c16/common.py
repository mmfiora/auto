# common.py
# Shared functions and constants for DBAASP API interaction

import csv
import requests
from config import Config

def load_ids(csv_file: str = None):
    """
    Read peptide IDs from CSV file.
    Accept a column named 'Peptide ID' or 'ID' (case-insensitive).
    Accept values like '51' or 'DBAASPS_51'.
    """
    if csv_file is None:
        csv_file = Config.INPUT_PEPTIDES_CSV
    
    with open(csv_file, encoding=Config.CSV_ENCODING) as f:
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
    resp = requests.get(Config.API_URL.format(id=pid), headers=Config.API_HEADERS, timeout=Config.API_TIMEOUT)
    resp.raise_for_status()
    return resp.json()