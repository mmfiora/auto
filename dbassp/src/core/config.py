# config.py
# Centralized configuration for DBAASP peptide analysis pipeline

import os
import logging
from dotenv import load_dotenv, find_dotenv

# Load .env from project root (safe no-op if file doesn't exist)
load_dotenv(find_dotenv(usecwd=True), override=False)

class Config:
    """Configuration settings with defaults and environment variable overrides."""
    
    # API Configuration
    API_URL = os.getenv("DBAASP_API_URL", "https://dbaasp.org/peptides/{id}")
    API_TIMEOUT = int(os.getenv("DBAASP_TIMEOUT", "20"))
    API_HEADERS = {
        "User-Agent": os.getenv("USER_AGENT", "Mozilla/5.0"),
        "Accept": "application/json"
    }
    
    # Symlinked folder paths (junctions → Google Drive)
    OUTPUT_DB_DIR     = os.getenv("OUTPUT_DB_DIR",     "output_db")
    PRESENTATIONS_DIR = os.getenv("PRESENTATIONS_DIR", "presentations")

    # File Paths (will be updated dynamically with Nterminus)
    INPUT_PEPTIDES_CSV    = os.getenv("INPUT_PEPTIDES_CSV",    "input/peptides.csv")
    OUTPUT_PHYSCHEM_CSV   = os.getenv("OUTPUT_PHYSCHEM_CSV",   f"{os.getenv('OUTPUT_DB_DIR', 'output_db')}/physchem.csv")
    OUTPUT_ACTIVITY_CSV   = os.getenv("OUTPUT_ACTIVITY_CSV",   f"{os.getenv('OUTPUT_DB_DIR', 'output_db')}/activity.csv")
    OUTPUT_NORMALIZED_CSV = os.getenv("OUTPUT_NORMALIZED_CSV", f"{os.getenv('OUTPUT_DB_DIR', 'output_db')}/activity_normalized.csv")
    OUTPUT_UNIFIED_CSV    = os.getenv("OUTPUT_UNIFIED_CSV",    f"{os.getenv('OUTPUT_DB_DIR', 'output_db')}/unified_results.csv")
    OUTPUT_LIPOPHILICITY_CSV = os.getenv("OUTPUT_LIPOPHILICITY_CSV", f"{os.getenv('OUTPUT_DB_DIR', 'output_db')}/lipophilicity.csv")
    OUTPUT_INTRINSIC_CSV  = os.getenv("OUTPUT_INTRINSIC_CSV",  f"{os.getenv('OUTPUT_DB_DIR', 'output_db')}/intrinsic_properties.csv")
    MIN_LIST_FILE         = os.getenv("MIN_LIST_FILE",         "input/list_min.txt")
    
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
                 logging.FileHandler("logs/pipeline.log", mode="a")
            ]
        )
        return logging.getLogger("dbaasp_pipeline")
    
    @classmethod
    def set_nterminus(cls, nterminus: str):
        """Set the Nterminus and update file paths accordingly."""
        cls._current_nterminus = nterminus
        d = cls.OUTPUT_DB_DIR
        cls.INPUT_PEPTIDES_CSV       = f"input/peptides_{nterminus}.csv"
        cls.OUTPUT_PHYSCHEM_CSV      = f"{d}/physchem_{nterminus}.csv"
        cls.OUTPUT_ACTIVITY_CSV      = f"{d}/activity_{nterminus}.csv"
        cls.OUTPUT_NORMALIZED_CSV    = f"{d}/activity_normalized_{nterminus}.csv"
        cls.OUTPUT_UNIFIED_CSV       = f"{d}/unified_results_{nterminus}.csv"
        cls.OUTPUT_LIPOPHILICITY_CSV = f"{d}/lipophilicity_{nterminus}.csv"
        cls.OUTPUT_INTRINSIC_CSV     = f"{d}/intrinsic_properties_{nterminus}.csv"
        cls.MIN_LIST_FILE            = f"input/list_min_{nterminus}.txt"
    
    @classmethod
    def get_nterminus(cls) -> str | None:
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