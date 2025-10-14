# physchem.py
import csv
import logging
from src.core import common
from src.core.config import Config
from src.core.exceptions import APIError, FileProcessingError

logger = logging.getLogger("dbaasp_pipeline")

def run():
    logger.info("Starting physicochemical data collection")
    
    try:
        ids = common.load_ids()
        if not ids:
            logger.warning(f"No IDs found in {Config.INPUT_PEPTIDES_CSV}")
            return

        logger.info(f"Processing {len(ids)} peptides for physicochemical data")
        data = []
        failed_peptides = []
        
        for i, pid in enumerate(ids, start=1):
            print(f"[{i}/{len(ids)}] Fetching physchem for peptide {pid} ...", end="", flush=True)
            try:
                d = common.fetch(pid)
                data.append(d)
                print(" ok")
            except APIError as e:
                logger.error(f"API error for peptide {pid}: {e}")
                failed_peptides.append(pid)
                print(f" fail (API error)")
            except Exception as e:
                logger.error(f"Unexpected error for peptide {pid}: {e}")
                failed_peptides.append(pid)
                print(f" fail (unexpected error)")

        if failed_peptides:
            logger.warning(f"Failed to fetch {len(failed_peptides)} peptides: {failed_peptides}")

        if not data:
            logger.error("No data fetched; nothing to write.")
            return

        logger.info(f"Successfully fetched data for {len(data)} peptides")
        
        # Collect property names preserving first-seen order, skip any property titled exactly "ID"
        props, seen = [], set()
        for d in data:
            for p in d.get("physicoChemicalProperties") or []:
                name = (p.get("name") or "").strip()
                if not name or name.upper() == "ID":
                    continue
                if name not in seen:
                    seen.add(name)
                    props.append(name)

        # Final header: base peptide columns + dynamic physchem properties
        header = ["Peptide ID", "N TERMINUS", "SEQUENCE", "C TERMINUS"] + props

        try:
            # Write CSV, always one row per peptide (fill missing properties with empty string)
            with open(Config.OUTPUT_PHYSCHEM_CSV, "w", newline="", encoding=Config.CSV_ENCODING) as f:
                w = csv.DictWriter(f, fieldnames=header)
                w.writeheader()
                for d in data:
                    row = {
                        "Peptide ID": str(d.get("id", "")),                          # numeric peptide ID
                        "N TERMINUS": (d.get("nTerminus") or {}).get("name", ""),    # only .name
                        "SEQUENCE": d.get("sequence", ""),                            # original casing
                        "C TERMINUS": (d.get("cTerminus") or {}).get("name", ""),    # only .name
                    }
                    # Map this peptide's properties
                    values = {}
                    for p in d.get("physicoChemicalProperties") or []:
                        name = (p.get("name") or "").strip()
                        if not name or name.upper() == "ID":
                            continue
                        value = str(p.get("value", "")).strip()
                        
                        # Adjust Net Charge for C-terminus peptides (C12, C16, etc.)
                        nterminus = Config.get_nterminus()
                        if name == "Net Charge" and nterminus and nterminus.startswith("C") and value:
                            try:
                                value = str(float(value) - 1.0)
                            except ValueError:
                                pass
                        
                        values[name] = value
                    # Fill all columns consistently
                    for name in props:
                        row[name] = values.get(name, "")
                    w.writerow(row)
                    
            logger.info(f"Successfully wrote physicochemical data to {Config.OUTPUT_PHYSCHEM_CSV}")
            
        except IOError as e:
            raise FileProcessingError(f"Error writing physicochemical file: {e}", filename=Config.OUTPUT_PHYSCHEM_CSV)
            
    except (FileProcessingError, APIError) as e:
        logger.error(f"Physicochemical collection failed: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in physicochemical collection: {e}")
        raise


