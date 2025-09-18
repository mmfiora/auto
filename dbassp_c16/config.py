# config.py
# Centralized configuration for DBAASP peptide analysis pipeline

import os

class Config:
    """Configuration settings with defaults and environment variable overrides."""
    
    # API Configuration
    API_URL = os.getenv("DBAASP_API_URL", "https://dbaasp.org/peptides/{id}")
    API_TIMEOUT = int(os.getenv("DBAASP_TIMEOUT", "20"))
    API_HEADERS = {
        "User-Agent": os.getenv("USER_AGENT", "Mozilla/5.0"),
        "Accept": "application/json"
    }
    
    # File Paths
    INPUT_PEPTIDES_CSV = os.getenv("INPUT_PEPTIDES_CSV", "peptides.csv")
    OUTPUT_PHYSCHEM_CSV = os.getenv("OUTPUT_PHYSCHEM_CSV", "physchem.csv")
    OUTPUT_ACTIVITY_CSV = os.getenv("OUTPUT_ACTIVITY_CSV", "activity.csv")
    OUTPUT_NORMALIZED_CSV = os.getenv("OUTPUT_NORMALIZED_CSV", "activity_normalized.csv")
    MIN_LIST_FILE = os.getenv("MIN_LIST_FILE", "list_min.txt")
    
    # File Encoding
    CSV_ENCODING = "utf-8-sig"
    
    # Molecular Weight Constants (Da)
    AA_MASS = {
        "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.15,
        "E": 129.12, "Q": 128.13, "G": 57.05,  "H": 137.14, "I": 113.16,
        "L": 113.16, "K": 128.17, "M": 131.20, "F": 147.18, "P": 97.12,
        "S": 87.08,  "T": 101.11, "W": 186.21, "Y": 163.18, "V": 99.13,
        "Z": 56.10,  # special C4 block
        "X": 110.0   # generic unknown residue
    }
    
    H2O_MASS = 18.02
    NTERM_MASS = {"C16": 239.2}     # N-terminus additions
    CTERM_MASS = {"AMD": -0.98}     # C-terminus modifications
    
    @classmethod
    def validate_files(cls) -> bool:
        """Check if required input files exist."""
        required_files = [cls.INPUT_PEPTIDES_CSV]
        missing = [f for f in required_files if not os.path.exists(f)]
        if missing:
            print(f"Missing required files: {missing}")
            return False
        return True