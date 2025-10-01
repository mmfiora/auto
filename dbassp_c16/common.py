# common.py
# Shared functions and constants for DBAASP API interaction

import csv
import glob
import logging
import requests
from config import Config
from exceptions import APIError, FileProcessingError, DataValidationError

logger = logging.getLogger("dbaasp_pipeline")

def auto_detect_nterminus() -> str:
    """
    Auto-detect Nterminus from available peptides_{Nterminus}.csv files.
    Returns the first Nterminus found, or raises an error if none found.
    """
    logger.info("Auto-detecting Nterminus from available peptides files")
    
    # Find all peptides_{Nterminus}.csv files
    peptide_files = glob.glob("peptides_*.csv")
    if not peptide_files:
        raise FileProcessingError(
            "No peptides_{Nterminus}.csv files found. Please ensure you have files like peptides_C16.csv",
            filename="peptides_*.csv"
        )
    
    # Sort files for consistent behavior
    peptide_files.sort()
    
    # Extract Nterminus from first filename (peptides_C16.csv -> C16)
    try:
        filename = peptide_files[0]
        nterminus = filename.split("_")[1].split(".")[0]
        
        if not nterminus:
            raise ValueError("Empty Nterminus")
            
        logger.info(f"Auto-detected Nterminus: {nterminus} from file: {filename}")
        
        # If multiple files exist, log them
        if len(peptide_files) > 1:
            logger.info(f"Found {len(peptide_files)} peptides files: {peptide_files}")
            logger.info(f"Using first file: {filename}")
            
        return nterminus
        
    except (IndexError, ValueError) as e:
        raise FileProcessingError(
            f"Invalid peptides filename format: {filename}. Expected format: peptides_{{Nterminus}}.csv",
            filename=filename
        )

def detect_nterminus_from_csv(csv_file: str) -> str:
    """
    Read the first row of peptides to detect the Nterminus value.
    Returns the Nterminus value found, or raises an error if none found.
    """
    logger.info(f"Detecting Nterminus from {csv_file}")
    
    try:
        with open(csv_file, encoding=Config.CSV_ENCODING) as f:
            r = csv.DictReader(f)
            
            # Check if N TERMINUS column exists
            if "N TERMINUS" not in (r.fieldnames or []):
                raise DataValidationError(
                    "No 'N TERMINUS' column found in CSV file",
                    field="N TERMINUS"
                )
            
            # Read first data row to get Nterminus
            for row in r:
                nterminus = (row.get("N TERMINUS") or "").strip()
                if nterminus:
                    logger.info(f"Detected Nterminus: {nterminus}")
                    return nterminus
                    
            raise DataValidationError(
                "No Nterminus value found in first row",
                field="N TERMINUS"
            )
            
    except FileNotFoundError:
        raise FileProcessingError(f"Input file not found: {csv_file}", filename=csv_file)
    except csv.Error as e:
        raise FileProcessingError(f"CSV parsing error: {e}", filename=csv_file)

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