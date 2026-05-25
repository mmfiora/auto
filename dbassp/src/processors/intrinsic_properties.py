# intrinsic_properties.py
# Combines physicochemical properties and min-list data into a single CSV
# representing the intrinsic molecular properties of peptides.
# Note: logP / logD (lipophilicity) are NOT included here — compute them
# in a separate project.

import csv
import re
import logging
from src.core.config import Config
from src.core.exceptions import FileProcessingError
from src.collectors.normalize_activity import load_min_map, get_z_prefix

logger = logging.getLogger("dbaasp_pipeline")


def create_intrinsic_csv(
    physchem_file: str | None = None,
    min_list_file: str | None = None,
    output_file: str | None = None
) -> None:
    """
    Create intrinsic properties CSV combining physicochemical data and min-list data.

    This output contains only the inherent molecular properties of peptides,
    independent of biological activity data. logP / logD are excluded
    (calculated externally).

    Args:
        physchem_file: Path to physchem.csv (default: Config.OUTPUT_PHYSCHEM_CSV)
        min_list_file: Path to list_min.txt (default: Config.MIN_LIST_FILE)
        output_file: Path for intrinsic properties output (default: intrinsic_properties.csv)
    """
    if physchem_file is None:
        physchem_file = Config.OUTPUT_PHYSCHEM_CSV
    if min_list_file is None:
        min_list_file = Config.MIN_LIST_FILE
    if output_file is None:
        output_file = Config.OUTPUT_INTRINSIC_CSV

    logger.info(f"Creating intrinsic properties CSV: {physchem_file} -> {output_file}")

    try:
        # Step 1: Load min list data using the shared space-delimited parser.
        # load_min_map() returns {NEW_SEQ -> (curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2)}
        min_list_map = load_min_map(min_list_file)
        MIN_LIST_COLS = ["curv_min", "npol_min", "ph_run", "npol_c0", "npol_c1", "npol_c2"]
        if min_list_map:
            min_list_headers = MIN_LIST_COLS
            logger.info(f"Loaded min list data for {len(min_list_map)} sequences")
        else:
            min_list_headers = []
            logger.warning(f"Min list file not found or empty: {min_list_file}")
            logger.warning("Continuing without min list data")

        # Step 2: Process physchem data and merge with min list
        intrinsic_rows = []

        with open(physchem_file, encoding=Config.CSV_ENCODING) as f:
            reader = csv.DictReader(f)

            if reader.fieldnames is None:
                raise FileProcessingError(f"No headers found in {physchem_file}", filename=physchem_file)

            # Create intrinsic header: physchem columns + derived scalar columns + min list columns
            physchem_headers = list(reader.fieldnames)
            new_cols = ["long_tail", "molecular_weight", "total_charge"]
            intrinsic_headers = physchem_headers + new_cols + min_list_headers

            row_count = 0
            min_matched = 0

            for row in reader:
                row_count += 1
                peptide_id = row["Peptide ID"]
                sequence   = (row.get("SEQUENCE") or "").upper()
                n_term     = row.get("N TERMINUS", "")
                c_term     = row.get("C TERMINUS", "")

                # Build the same NEW_SEQ key used by normalize_activity
                z_prefix = get_z_prefix(n_term)
                new_seq  = f"{z_prefix}{sequence}{'01' if (c_term or '').upper() == 'AMD' else '00'}"

                # ── Derived scalars ──────────────────────────────────────────────────

                # 1. long_tail: chain length extracted from N TERMINUS (e.g. C12 -> "12")
                long_tail = ""
                if n_term:
                    match = re.search(r"C(\d+)", n_term.upper())
                    if match:
                        long_tail = match.group(1)

                # 2. molecular_weight: sequence residues + N-term fatty acid + C-term mod
                mw = sum(Config.AA_MASS.get(a.upper(), 110.0) for a in (sequence or "")) + Config.H2O_MASS
                if n_term and n_term.upper() in Config.NTERM_MASS:
                    mw += Config.NTERM_MASS[n_term.upper()]
                if c_term and c_term.upper() in Config.CTERM_MASS:
                    mw += Config.CTERM_MASS[c_term.upper()]

                # 3. total_charge: Net Charge from DBAASP adjusted for acylation / amidation
                #    Acylated N-term (C4–C20) loses +1; amidated C-term (AMD) gains +1
                raw_net_charge = row.get("Net Charge", "")
                total_charge = ""
                if raw_net_charge != "":
                    try:
                        tc = float(raw_net_charge)
                        if long_tail:          # fatty-acid acylation removes free amine
                            tc -= 1.0
                        if (c_term or "").upper() == "AMD":
                            tc += 1.0
                        total_charge = f"{tc:.2f}"
                    except ValueError:
                        total_charge = ""

                # Build row
                intrinsic_row = dict(row)
                intrinsic_row["long_tail"] = long_tail
                intrinsic_row["molecular_weight"] = f"{mw:.2f}"
                intrinsic_row["total_charge"] = total_charge

                # Add min list properties — exact NEW_SEQ match only.
                # 00 (AMD) and 01 (non-AMD) are different simulations; never cross-assign.
                if new_seq in min_list_map:
                    curv, npol, ph, nc0, nc1, nc2 = min_list_map[new_seq]
                    intrinsic_row["curv_min"] = curv
                    intrinsic_row["npol_min"] = npol
                    intrinsic_row["ph_run"]   = ph
                    intrinsic_row["npol_c0"]  = nc0
                    intrinsic_row["npol_c1"]  = nc1
                    intrinsic_row["npol_c2"]  = nc2
                    min_matched += 1
                else:
                    for header in min_list_headers:
                        intrinsic_row[header] = ""

                intrinsic_rows.append(intrinsic_row)

        logger.info(f"Processed {row_count} peptides")
        logger.info(f"Matched {min_matched} with min list data")

        # Step 3: Write intrinsic properties CSV
        with open(output_file, "w", encoding=Config.CSV_ENCODING, newline="") as f:
            writer = csv.DictWriter(f, fieldnames=intrinsic_headers)
            writer.writeheader()
            writer.writerows(intrinsic_rows)

        logger.info(f"Successfully created intrinsic properties CSV: {output_file}")
        logger.info(f"Final CSV: {len(intrinsic_rows)} rows × {len(intrinsic_headers)} columns")
        logger.info("Intrinsic properties CSV includes:")
        logger.info("  - Basic peptide info (ID, N-terminus, sequence, C-terminus)")
        logger.info(f"  - Physicochemical properties from DBAASP ({len(physchem_headers) - 4} columns)")
        logger.info("  - Derived: long_tail, molecular_weight, total_charge")
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
