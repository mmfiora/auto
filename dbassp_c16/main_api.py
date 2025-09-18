# main_api.py
import time
import logging
from config import Config
from exceptions import DBAASSPError
import physchem
import activity
import normalize_activity
import unified_results

def main():
    # Setup logging
    logger = Config.setup_logging()
    logger.info("Starting DBAASP pipeline")
    
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

