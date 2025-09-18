# common.py
# Shared functions and constants for DBAASP API interaction

import csv
import logging
import requests
from config import Config
from exceptions import APIError, FileProcessingError, DataValidationError

logger = logging.getLogger("dbaasp_pipeline")

def load_ids(csv_file: str | None = None):
    """
    Read peptide IDs from CSV file.
    Accept a column named 'Peptide ID' or 'ID' (case-insensitive).
    Accept values like '51' or 'DBAASPS_51'.
    """
    if csv_file is None:
        csv_file = Config.INPUT_PEPTIDES_CSV
    
    logger.info(f"Loading peptide IDs from {csv_file}")
    
    try:
        with open(csv_file, encoding=Config.CSV_ENCODING) as f:
            r = csv.DictReader(f)
            col = None
            for h in r.fieldnames or []:
                if h.lower() in ("peptide id", "id"):
                    col = h
                    break
            if not col:
                raise DataValidationError(
                    "No valid ID column found. Expected 'Peptide ID' or 'ID'",
                    field="column_names"
                )
            
            ids = []
            for row_num, row in enumerate(r, start=2):  # Start at 2 for header
                raw = (row.get(col) or "").strip()
                if not raw:
                    logger.warning(f"Empty ID in row {row_num}, skipping")
                    continue
                try:
                    ids.append(int(raw.split("_")[-1]))
                except ValueError as e:
                    logger.warning(f"Invalid ID format '{raw}' in row {row_num}: {e}")
                    continue
                    
            logger.info(f"Successfully loaded {len(ids)} peptide IDs")
            return ids
            
    except FileNotFoundError:
        raise FileProcessingError(f"Input file not found: {csv_file}", filename=csv_file)
    except csv.Error as e:
        raise FileProcessingError(f"CSV parsing error: {e}", filename=csv_file)

def fetch(pid: int):
    """Fetch a single peptide JSON from DBAASP API."""
    logger.debug(f"Fetching peptide {pid} from API")
    
    try:
        resp = requests.get(
            Config.API_URL.format(id=pid), 
            headers=Config.API_HEADERS, 
            timeout=Config.API_TIMEOUT
        )
        resp.raise_for_status()
        data = resp.json()
        logger.debug(f"Successfully fetched peptide {pid}")
        return data
        
    except requests.exceptions.Timeout:
        raise APIError(f"Timeout fetching peptide {pid}", peptide_id=pid)
    except requests.exceptions.ConnectionError:
        raise APIError(f"Connection error fetching peptide {pid}", peptide_id=pid)
    except requests.exceptions.HTTPError as e:
        raise APIError(
            f"HTTP {e.response.status_code} error fetching peptide {pid}: {e}", 
            peptide_id=pid, 
            status_code=e.response.status_code
        )
    except requests.exceptions.JSONDecodeError:
        raise APIError(f"Invalid JSON response for peptide {pid}", peptide_id=pid)