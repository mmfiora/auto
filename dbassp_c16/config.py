# config.py
# Centralized configuration for DBAASP peptide analysis pipeline

import os
import logging

class Config:
    """Configuration settings with defaults and environment variable overrides."""
    
    # API Configuration
    API_URL = os.getenv("DBAASP_API_URL", "https://dbaasp.org/peptides/{id}")
    API_TIMEOUT = int(os.getenv("DBAASP_TIMEOUT", "20"))
    API_HEADERS = {
        "User-Agent": os.getenv("USER_AGENT", "Mozilla/5.0"),
        "Accept": "application/json"
    }
    
    # File Paths (will be updated dynamically with Nterminus)
    INPUT_PEPTIDES_CSV = os.getenv("INPUT_PEPTIDES_CSV", "peptides.csv")
    OUTPUT_PHYSCHEM_CSV = os.getenv("OUTPUT_PHYSCHEM_CSV", "physchem.csv")
    OUTPUT_ACTIVITY_CSV = os.getenv("OUTPUT_ACTIVITY_CSV", "activity.csv")
    OUTPUT_NORMALIZED_CSV = os.getenv("OUTPUT_NORMALIZED_CSV", "activity_normalized.csv")
    OUTPUT_UNIFIED_CSV = os.getenv("OUTPUT_UNIFIED_CSV", "unified_results.csv")
    MIN_LIST_FILE = os.getenv("MIN_LIST_FILE", "list_min.txt")
    
    # Current Nterminus (set dynamically during pipeline execution)
    _current_nterminus = None
    
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
    
    # Logging Configuration
    LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO")
    LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    
    @classmethod
    def setup_logging(cls):
        """Setup centralized logging configuration."""
        logging.basicConfig(
            level=getattr(logging, cls.LOG_LEVEL.upper()),
            format=cls.LOG_FORMAT,
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler("pipeline.log", mode="a")
            ]
        )
        return logging.getLogger("dbaasp_pipeline")
    
    @classmethod
    def set_nterminus(cls, nterminus: str):
        """Set the Nterminus and update file paths accordingly."""
        cls._current_nterminus = nterminus
        cls.INPUT_PEPTIDES_CSV = f"peptides_{nterminus}.csv"
        cls.OUTPUT_PHYSCHEM_CSV = f"physchem_{nterminus}.csv"
        cls.OUTPUT_ACTIVITY_CSV = f"activity_{nterminus}.csv"
        cls.OUTPUT_NORMALIZED_CSV = f"activity_normalized_{nterminus}.csv"
        cls.OUTPUT_UNIFIED_CSV = f"unified_results_{nterminus}.csv"
        cls.MIN_LIST_FILE = f"list_min_{nterminus}.txt"
    
    @classmethod
    def get_nterminus(cls) -> str:
        """Get the current Nterminus."""
        return cls._current_nterminus
    
    @classmethod
    def validate_files(cls) -> bool:
        """Check if required input files exist."""
        logger = logging.getLogger("dbaasp_pipeline")
        required_files = [cls.INPUT_PEPTIDES_CSV]
        missing = [f for f in required_files if not os.path.exists(f)]
        if missing:
            logger.error(f"Missing required files: {missing}")
            return False
        logger.info(f"All required files found: {required_files}")
        return True