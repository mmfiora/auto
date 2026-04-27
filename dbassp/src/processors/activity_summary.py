# activity_summary.py
# Produces a slim summary CSV from activity_normalized with only the
# columns needed for downstream analysis.

import csv
import logging
from src.core.config import Config
from src.core.exceptions import FileProcessingError

logger = logging.getLogger("dbaasp_pipeline")

SUMMARY_COLS = [
    "N TERMINUS",
    "SEQUENCE",
    "C TERMINUS",
    "lower_ugml",
    "upper_ugml",
    "species",
    "strain",
    "medium",
    "curv_min",
    "npol_min",
    "ph_run",
    "activity",
]


def run(infile: str | None = None, outfile: str | None = None) -> None:
    """
    Read activity_normalized CSV and write a slim summary with only the
    columns: N TERMINUS, SEQUENCE, C TERMINUS, lower_ugml, upper_ugml,
    species, strain, medium, curv_min, npol_min, ph_run, activity.
    """
    if infile is None:
        infile = Config.OUTPUT_NORMALIZED_CSV
    if outfile is None:
        nt = Config.get_nterminus() or "unknown"
        outfile = f"data/output/activity_summary_{nt}.csv"

    logger.info(f"Creating activity summary CSV: {infile} -> {outfile}")

    try:
        with open(infile, encoding=Config.CSV_ENCODING) as f:
            reader = csv.DictReader(f)

            # Warn about any requested columns that are missing
            missing = [c for c in SUMMARY_COLS if c not in (reader.fieldnames or [])]
            if missing:
                logger.warning(f"Columns not found in input and will be blank: {missing}")

            rows = []
            skipped = 0
            for row in reader:
                if not row.get("curv_min", ""):
                    skipped += 1
                    continue
                rows.append({col: row.get(col, "") for col in SUMMARY_COLS})

        with open(outfile, "w", encoding=Config.CSV_ENCODING, newline="") as f:
            writer = csv.DictWriter(f, fieldnames=SUMMARY_COLS)
            writer.writeheader()
            writer.writerows(rows)

        logger.info(f"Successfully wrote activity summary: {outfile} ({len(rows)} rows, {skipped} skipped — no curv_min)")

    except FileNotFoundError:
        raise FileProcessingError(f"Input file not found: {infile}", filename=infile)
    except csv.Error as e:
        raise FileProcessingError(f"CSV processing error: {e}", filename=infile)
    except IOError as e:
        raise FileProcessingError(f"File I/O error: {e}", filename=outfile)
    except Exception as e:
        logger.error(f"Unexpected error creating activity summary: {e}")
        raise
