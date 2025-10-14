# main_api.py
import argparse
import os
import time
import logging
import glob
import re
from src.core.config import Config
from src.core.exceptions import DBAASSPError
from src.collectors import physchem
from src.collectors import activity
from src.collectors import normalize_activity
from src.processors import unified_results
from src.core import common
import os
os.makedirs("logs", exist_ok=True)
def main():
    parser = argparse.ArgumentParser(
        description="DBAASP peptide analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 main_api.py                    # Interactive menu
  python3 main_api.py --nterminus C16    # Use specific Nterminus
        """
    )
    parser.add_argument(
        '--nterminus', 
        type=str, 
        help='Specify Nterminus to use (e.g., C16, C12). If not specified, shows interactive menu'
    )
    args = parser.parse_args()
    
    logger = Config.setup_logging()
    logger.info("Starting DBAASP pipeline")
    
    try:
        if args.nterminus:
            nterminus = args.nterminus
            expected_file = f"data/input/peptides_{nterminus}.csv"
            
            if not os.path.exists(expected_file):
                logger.error(f"Specified peptides file not found: {expected_file}")
                print(f"Error: File {expected_file} not found")
                return 1
                
            Config.set_nterminus(nterminus)
            logger.info(f"Using specified Nterminus: {nterminus}")
            print(f"Using Nterminus: {nterminus}")
            
        else:
            peptide_files = glob.glob("data/input/peptides_*.csv")
            if not peptide_files:
                logger.error("No peptides_*.csv files found in data/input/")
                print("Error: No peptides files found in data/input/")
                return 1
            
            nterminus_list = []
            pattern = re.compile(r'peptides_([A-Z]+\d+)\.csv$')
            
            for f in sorted(peptide_files):
                basename = os.path.basename(f)
                match = pattern.match(basename)
                if match:
                    nt = match.group(1)
                    if nt not in nterminus_list:
                        nterminus_list.append(nt)
            
            if not nterminus_list:
                print("Error: No valid peptides files found")
                return 1
            
            print("\nAvailable Nterminus:")
            for i, nt in enumerate(nterminus_list, 1):
                print(f"{i}. {nt}")
            
            choice = input("\nSelect Nterminus number: ").strip()
            if not choice.isdigit() or int(choice) < 1 or int(choice) > len(nterminus_list):
                print("Invalid selection")
                return 1
            
            nterminus = nterminus_list[int(choice) - 1]
            Config.set_nterminus(nterminus)
            logger.info(f"Selected Nterminus: {nterminus}")
            print(f"\nUsing Nterminus: {nterminus}")
        
        print(f"Input: {Config.INPUT_PEPTIDES_CSV}")
        print(f"Output: {Config.OUTPUT_PHYSCHEM_CSV}, {Config.OUTPUT_ACTIVITY_CSV}, {Config.OUTPUT_UNIFIED_CSV}\n")
        
    except Exception as e:
        logger.error(f"Failed to configure Nterminus: {e}")
        print(f"Error: {e}")
        return 1
    
    if not Config.validate_files():
        logger.error("Configuration validation failed")
        return 1
    
    t0 = time.perf_counter()

    try:
        logger.info("Running physchem collection...")
        print("Running physchem ...")
        physchem.run()

        logger.info("Running activity collection...")
        print("Running activity ...")
        activity.run()

        logger.info("Normalizing activity data...")
        print("Normalizing activity ...")
        normalize_activity.run()

        logger.info("Creating unified results CSV...")
        print("Creating unified results ...")
        unified_results.run()

        t1 = time.perf_counter()
        elapsed = t1 - t0
        logger.info(f"Pipeline completed successfully in {elapsed:.2f} seconds")
        print(f"\nDone! Total time: {elapsed:.2f} s")
        return 0
        
    except DBAASSPError as e:
        logger.error(f"Pipeline failed with DBAASP error: {e}")
        print(f"Error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Pipeline failed with unexpected error: {e}")
        print(f"Unexpected error: {e}")
        return 1

if __name__ == "__main__":
    exit(main())
