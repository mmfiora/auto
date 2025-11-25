# intrinsic_properties.py
# Combines physicochemical and lipophilicity properties into a single CSV
# representing the intrinsic molecular properties of peptides

import csv
import logging
from src.core.config import Config
from src.core.exceptions import FileProcessingError

logger = logging.getLogger("dbaasp_pipeline")

def create_intrinsic_csv(
    physchem_file: str | None = None,
    lipophilicity_file: str | None = None,
    min_list_file: str | None = None,
    output_file: str | None = None
) -> None:
    """
    Create intrinsic properties CSV combining physicochemical and lipophilicity data.
    
    This output contains only the inherent molecular properties of peptides,
    independent of biological activity data.
    
    Args:
        physchem_file: Path to physchem.csv (default: Config.OUTPUT_PHYSCHEM_CSV)
        lipophilicity_file: Path to lipophilicity.csv (default: Config.OUTPUT_LIPOPHILICITY_CSV)
        min_list_file: Path to list_min.txt (default: Config.MIN_LIST_FILE)
        output_file: Path for intrinsic properties output (default: intrinsic_properties.csv)
    """
    if physchem_file is None:
        physchem_file = Config.OUTPUT_PHYSCHEM_CSV
    if lipophilicity_file is None:
        lipophilicity_file = Config.OUTPUT_LIPOPHILICITY_CSV
    if min_list_file is None:
        min_list_file = Config.MIN_LIST_FILE
    if output_file is None:
        output_file = Config.OUTPUT_INTRINSIC_CSV
    
    logger.info(f"Creating intrinsic properties CSV: {physchem_file} + {lipophilicity_file} -> {output_file}")
    
    try:
        # Step 1: Load lipophilicity data by Peptide ID
        lipophilicity_data = {}
        lipophilicity_headers = []
        
        try:
            with open(lipophilicity_file, encoding=Config.CSV_ENCODING) as f:
                reader = csv.DictReader(f)
                if reader.fieldnames is not None:
                    # Exclude basic peptide info, keep only lipophilicity properties
                    lipophilicity_headers = [h for h in reader.fieldnames 
                                            if h not in ["Peptide ID", "N TERMINUS", "SEQUENCE", "SMILES", "C TERMINUS"]]
                    
                    for row in reader:
                        peptide_id = row["Peptide ID"]
                        lipophilicity_props = {h: row.get(h, "") for h in lipophilicity_headers}
                        lipophilicity_data[peptide_id] = lipophilicity_props
            
            logger.info(f"Loaded lipophilicity data for {len(lipophilicity_data)} peptides")
            logger.info(f"Lipophilicity properties: {lipophilicity_headers}")
        except FileNotFoundError:
            logger.warning(f"Lipophilicity file not found: {lipophilicity_file}")
            logger.warning("Continuing without lipophilicity data")
            lipophilicity_headers = []
        
        # Step 2: Load min list data if available
        min_list_data = {}
        min_list_headers = []
        
        try:
            with open(min_list_file, encoding=Config.CSV_ENCODING) as f:
                # Skip first line if it's a header or read it
                lines = f.readlines()
                if lines:
                    # Check if first line contains headers
                    first_line = lines[0].strip()
                    if 'sequence' in first_line.lower() or 'pH' in first_line:
                        # Parse as CSV with headers
                        reader = csv.DictReader(lines)
                        if reader.fieldnames is not None:
                            # Exclude sequence column, keep min list properties
                            min_list_headers = [h for h in reader.fieldnames 
                                              if h.lower() not in ["sequence"]]
                            
                            for row in reader:
                                sequence = row.get("sequence", "").replace("ZZZ", "").replace("00", "").replace("01", "")
                                min_list_props = {h: row.get(h, "") for h in min_list_headers}
                                min_list_data[sequence] = min_list_props
                        
                        logger.info(f"Loaded min list data for {len(min_list_data)} sequences")
                        logger.info(f"Min list properties: {min_list_headers}")
        except FileNotFoundError:
            logger.warning(f"Min list file not found: {min_list_file}")
            logger.warning("Continuing without min list data")
            min_list_headers = []
        
        # Step 3: Process physchem data and merge with lipophilicity and min list
        intrinsic_rows = []
        
        with open(physchem_file, encoding=Config.CSV_ENCODING) as f:
            reader = csv.DictReader(f)
            
            if reader.fieldnames is None:
                raise FileProcessingError(f"No headers found in {physchem_file}", filename=physchem_file)
            
            # Create intrinsic header: physchem columns + lipophilicity columns + min list columns
            physchem_headers = list(reader.fieldnames)
            intrinsic_headers = physchem_headers + lipophilicity_headers + min_list_headers
            
            row_count = 0
            lipo_matched = 0
            min_matched = 0
            
            for row in reader:
                row_count += 1
                peptide_id = row["Peptide ID"]
                sequence = row.get("SEQUENCE", "")
                
                # Create intrinsic row starting with physchem data
                intrinsic_row = dict(row)
                
                # Add lipophilicity properties if available
                if peptide_id in lipophilicity_data:
                    intrinsic_row.update(lipophilicity_data[peptide_id])
                    lipo_matched += 1
                else:
                    # Fill missing lipophilicity data with empty strings
                    for header in lipophilicity_headers:
                        intrinsic_row[header] = ""
                
                # Add min list properties if available (match by sequence)
                if sequence in min_list_data:
                    intrinsic_row.update(min_list_data[sequence])
                    min_matched += 1
                else:
                    # Fill missing min list data with empty strings
                    for header in min_list_headers:
                        intrinsic_row[header] = ""
                
                intrinsic_rows.append(intrinsic_row)
        
        logger.info(f"Processed {row_count} peptides")
        logger.info(f"Matched {lipo_matched} with lipophilicity data")
        logger.info(f"Matched {min_matched} with min list data")
        
        # Step 4: Write intrinsic properties CSV
        with open(output_file, "w", encoding=Config.CSV_ENCODING, newline="") as f:
            writer = csv.DictWriter(f, fieldnames=intrinsic_headers)
            writer.writeheader()
            writer.writerows(intrinsic_rows)
        
        logger.info(f"Successfully created intrinsic properties CSV: {output_file}")
        logger.info(f"Final CSV contains {len(intrinsic_rows)} rows and {len(intrinsic_headers)} columns")
        
        # Log summary of what's included
        logger.info("Intrinsic properties CSV includes:")
        logger.info("  - Basic peptide info (ID, N-terminus, sequence, C-terminus)")
        logger.info(f"  - Physicochemical properties ({len(physchem_headers) - 4} columns)")
        logger.info(f"  - Lipophilicity properties ({len(lipophilicity_headers)} columns)")
        logger.info(f"  - Min list properties ({len(min_list_headers)} columns)")
        
    except FileNotFoundError as e:
        raise FileProcessingError(f"Input file not found: {e}", filename=str(e))
    except csv.Error as e:
        raise FileProcessingError(f"CSV processing error: {e}")
    except IOError as e:
        raise FileProcessingError(f"File I/O error: {e}", filename=output_file)
    except Exception as e:
        logger.error(f"Unexpected error creating intrinsic properties CSV: {e}")
        raise

def run(output_file: str | None = None) -> None:
    """Convenience function to create intrinsic properties CSV with default parameters."""
    create_intrinsic_csv(output_file=output_file)

if __name__ == "__main__":
    run()
