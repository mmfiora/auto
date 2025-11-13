# activity.py
import csv
import logging
from src.core import common
from src.core.config import Config
from src.core.exceptions import APIError, FileProcessingError

logger = logging.getLogger("dbaasp_pipeline")

def get_activities(peptide_json):
    a = peptide_json.get("targetActivities")
    if a is None:
        a = peptide_json.get("activityAgainstTargetSpecies")
    return a if isinstance(a, list) else []

def get_unusual_amino_acids(peptide_json):
    unusual = peptide_json.get("unusualAminoAcids", [])
    if not unusual:
        return ""
    names = [aa.get("modificationType", {}).get("name", "") for aa in unusual if aa.get("modificationType")]
    return ", ".join(filter(None, names))

def collect_activity_keys(all_peptides):
    exclude = {"id", "activity", "activityMeasureValue"}
    order, seen = [], set()
    for d in all_peptides:
        for a in get_activities(d):
            for k in a.keys():
                if not k or k.lower() in exclude:
                    continue
                if k not in seen:
                    seen.add(k)
                    order.append(k)
    return order

def flatten_value(v):
    if v is None:
        return ""
    if isinstance(v, dict):
        return v.get("name", str(v))
    return v

def run():
    logger.info("Starting activity data collection")
    
    try:
        ids = common.load_ids()
        if not ids:
            logger.warning(f"No IDs found in {Config.INPUT_PEPTIDES_CSV}")
            return

        logger.info(f"Processing {len(ids)} peptides for activity data")
        peptides = []
        failed_peptides = []
        
        for i, pid in enumerate(ids, start=1):
            print(f"[{i}/{len(ids)}] Fetching activity for peptide {pid} ...", end="", flush=True)
            try:
                peptides.append(common.fetch(pid))
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

        if not peptides:
            logger.error("No data fetched; nothing to write.")
            return

        logger.info(f"Successfully fetched data for {len(peptides)} peptides")
        activity_cols = collect_activity_keys(peptides)
        header = ["Peptide ID", "N TERMINUS", "SEQUENCE", "C TERMINUS", "Unusual Amino Acids"] + activity_cols

        try:
            with open(Config.OUTPUT_ACTIVITY_CSV, "w", newline="", encoding=Config.CSV_ENCODING) as f:
                w = csv.DictWriter(f, fieldnames=header)
                w.writeheader()

                for d in peptides:
                    base = {
                        "Peptide ID": str(d.get("id", "")),
                        "N TERMINUS": (d.get("nTerminus") or {}).get("name", ""),
                        "SEQUENCE": d.get("sequence", ""),
                        "C TERMINUS": (d.get("cTerminus") or {}).get("name", ""),
                        "Unusual Amino Acids": get_unusual_amino_acids(d),
                    }
                    acts = get_activities(d)
                    if not acts:
                        w.writerow(base)
                        continue
                    for a in acts:
                        row = dict(base)
                        for k in activity_cols:
                            row[k] = flatten_value(a.get(k))
                        w.writerow(row)
                        
            logger.info(f"Successfully wrote activity data to {Config.OUTPUT_ACTIVITY_CSV}")
            
        except IOError as e:
            raise FileProcessingError(f"Error writing activity file: {e}", filename=Config.OUTPUT_ACTIVITY_CSV)
            
    except (FileProcessingError, APIError) as e:
        logger.error(f"Activity collection failed: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in activity collection: {e}")
        raise
