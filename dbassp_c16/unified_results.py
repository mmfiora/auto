# unified_results.py
# Combines physchem, activity, and normalized data into a single CSV

import csv
import logging
from config import Config
from exceptions import FileProcessingError

logger = logging.getLogger("dbaasp_pipeline")

def create_unified_csv(
    physchem_file: str | None = None,
    normalized_file: str | None = None,
    output_file: str | None = None
) -> None:
    """
    Create unified CSV combining physicochemical properties with normalized activity data.
    
    Args:
        physchem_file: Path to physchem.csv (default: Config.OUTPUT_PHYSCHEM_CSV)
        normalized_file: Path to activity_normalized.csv (default: Config.OUTPUT_NORMALIZED_CSV)
        output_file: Path for unified output (default: unified_results.csv)
    """
    if physchem_file is None:
        physchem_file = Config.OUTPUT_PHYSCHEM_CSV
    if normalized_file is None:
        normalized_file = Config.OUTPUT_NORMALIZED_CSV
    if output_file is None:
        output_file = Config.OUTPUT_UNIFIED_CSV
    
    logger.info(f"Creating unified CSV: {physchem_file} + {normalized_file} -> {output_file}")
    
    try:
        # Step 1: Load physicochemical properties by Peptide ID
        physchem_data = {}
        physchem_headers = []
        
        with open(physchem_file, encoding=Config.CSV_ENCODING) as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise FileProcessingError(f"No headers found in {physchem_file}", filename=physchem_file)
            physchem_headers = [h for h in reader.fieldnames if h not in ["Peptide ID", "N TERMINUS", "SEQUENCE", "C TERMINUS"]]
            
            for row in reader:
                peptide_id = row["Peptide ID"]
                # Store only the physchem-specific columns (exclude basic peptide info)
                physchem_props = {h: row.get(h, "") for h in physchem_headers}
                physchem_data[peptide_id] = physchem_props
        
        logger.info(f"Loaded physicochemical data for {len(physchem_data)} peptides")
        logger.info(f"Physicochemical properties: {len(physchem_headers)} columns")
        
        # Step 2: Process normalized activity data and merge with physchem
        unified_rows = []
        
        with open(normalized_file, encoding=Config.CSV_ENCODING) as f:
            reader = csv.DictReader(f)
            
            # Create unified header: normalized activity columns + physchem columns
            if reader.fieldnames is None:
                raise FileProcessingError(f"No headers found in {normalized_file}", filename=normalized_file)
            activity_headers = list(reader.fieldnames)
            unified_headers = activity_headers + physchem_headers
            
            row_count = 0
            matched_count = 0
            
            for row in reader:
                row_count += 1
                peptide_id = row["Peptide ID"]
                
                # Create unified row starting with activity/normalized data
                unified_row = dict(row)
                
                # Add physicochemical properties if available
                if peptide_id in physchem_data:
                    unified_row.update(physchem_data[peptide_id])
                    matched_count += 1
                else:
                    # Fill missing physchem data with empty strings
                    for header in physchem_headers:
                        unified_row[header] = ""
                    logger.warning(f"No physicochemical data found for Peptide ID {peptide_id}")
                
                unified_rows.append(unified_row)
        
        logger.info(f"Processed {row_count} activity rows, matched {matched_count} with physchem data")
        
        # Step 3: Write unified CSV
        with open(output_file, "w", encoding=Config.CSV_ENCODING, newline="") as f:
            writer = csv.DictWriter(f, fieldnames=unified_headers)
            writer.writeheader()
            writer.writerows(unified_rows)
        
        logger.info(f"Successfully created unified CSV: {output_file}")
        logger.info(f"Final CSV contains {len(unified_rows)} rows and {len(unified_headers)} columns")
        
        # Log summary of what's included
        logger.info("Unified CSV includes:")
        logger.info("  - Basic peptide info (ID, sequence, terminus)")
        logger.info("  - Activity data (target species, MIC, concentration)")
        logger.info("  - Normalized concentrations (µM conversion)")
        logger.info("  - Min list data (curv_min, npol_min)")
        logger.info(f"  - Physicochemical properties ({len(physchem_headers)} columns)")
        
    except FileNotFoundError as e:
        raise FileProcessingError(f"Input file not found: {e}", filename=str(e))
    except csv.Error as e:
        raise FileProcessingError(f"CSV processing error: {e}")
    except IOError as e:
        raise FileProcessingError(f"File I/O error: {e}", filename=output_file)
    except Exception as e:
        logger.error(f"Unexpected error creating unified CSV: {e}")
        raise

def run(output_file: str | None = None) -> None:
    """Convenience function to create unified CSV with default parameters."""
    create_unified_csv(output_file=output_file)

if __name__ == "__main__":
    run()