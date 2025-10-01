# main_api.py
import argparse
import os
import time
import logging
from config import Config
from exceptions import DBAASSPError
import physchem
import activity
import normalize_activity
import unified_results
import common

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="DBAASP peptide analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 main_api.py                    # Auto-detect Nterminus from available files
  python3 main_api.py --nterminus C16    # Use specific Nterminus (peptides_C16.csv)
  python3 main_api.py --nterminus C12    # Use specific Nterminus (peptides_C12.csv)
        """
    )
    parser.add_argument(
        '--nterminus', 
        type=str, 
        help='Specify Nterminus to use (e.g., C16, C12). If not specified, auto-detects from available peptides_*.csv files'
    )
    args = parser.parse_args()
    
    # Setup logging
    logger = Config.setup_logging()
    logger.info("Starting DBAASP pipeline")
    
    # Determine Nterminus: specified or auto-detected
    try:
        if args.nterminus:
            # User specified Nterminus
            nterminus = args.nterminus
            expected_file = f"peptides_{nterminus}.csv"
            
            if not os.path.exists(expected_file):
                logger.error(f"Specified peptides file not found: {expected_file}")
                print(f"Error: File {expected_file} not found")
                print(f"Please ensure {expected_file} exists or use auto-detection (run without --nterminus)")
                return 1
                
            Config.set_nterminus(nterminus)
            logger.info(f"Using specified Nterminus: {nterminus}")
            print(f"Using specified Nterminus: {nterminus}")
            
        else:
            # Auto-detect Nterminus
            nterminus = common.auto_detect_nterminus()
            Config.set_nterminus(nterminus)
            logger.info(f"Auto-detected Nterminus: {nterminus}")
            print(f"Auto-detected Nterminus: {nterminus}")
        
        # Show configured files
        print(f"Using input files: {Config.INPUT_PEPTIDES_CSV}, {Config.MIN_LIST_FILE}")
        print(f"Output files will be: {Config.OUTPUT_PHYSCHEM_CSV}, {Config.OUTPUT_ACTIVITY_CSV}, {Config.OUTPUT_NORMALIZED_CSV}")
        
    except Exception as e:
        logger.error(f"Failed to configure Nterminus: {e}")
        print(f"Error configuring Nterminus: {e}")
        if not args.nterminus:
            print("Please ensure you have peptides_{Nterminus}.csv files available (e.g., peptides_C16.csv)")
        return 1
    
    # Validate configuration
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
        print(f"Done. Total elapsed: {elapsed:.2f} s")
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

