import glob
import re
import os
import sys

# Add current dir to path to import src
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from src.core.config import Config
from src.collectors import normalize_activity
from src.processors import activity_summary
from src.processors import unified_results
from src.processors import intrinsic_properties
import logging

logger = Config.setup_logging()

peptide_files = glob.glob("data/input/peptides_*.csv")
pattern = re.compile(r'peptides_([A-Za-z]*\d+)\.csv$', re.IGNORECASE)

for f in peptide_files:
    basename = os.path.basename(f)
    match = pattern.match(basename)
    if match:
        nt = match.group(1)
        print(f"--- Processing N-terminus: {nt} ---")
        Config.set_nterminus(nt)
        if os.path.exists(Config.OUTPUT_ACTIVITY_CSV):
            normalize_activity.run()
            activity_summary.run()
            unified_results.run()
            intrinsic_properties.run()
        else:
            print(f"Skipping {nt}, no {Config.OUTPUT_ACTIVITY_CSV} found")
