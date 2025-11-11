# normalize_activity.py
# Minimal normalizer:
# - Compute peptide MW from SEQUENCE (+ optional N- and C-terminus mods)
# - Split concentration into lower/upper (handles "a-b", ">x", "<x", single value)
# - Convert each bound to µM
# - Split targetSpecies into species and strain using a simple heuristic
# - Write activity_normalized.csv (or a custom outfile)
# + Join por NEW_SEQ con list_min.txt para curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2

import csv
import logging
import re
from src.core.config import Config
from src.core.exceptions import FileProcessingError, DataValidationError

logger = logging.getLogger("dbaasp_pipeline")

# MIC Threshold Configuration
# Sets the boundary for classifying peptides as having significant MIC activity
# MIC values <= this threshold are classified as 1 (has MIC activity)
# MIC values > this threshold or missing are classified as 0 (no MIC activity)
MIC_THRESHOLD_uM = 25.0  # microMolar threshold for MIC classification

def classify_mic_activity(lower_uM: str, upper_uM: str) -> int:
    """
    Classify MIC activity based on upper bound concentration value.
    
    Returns 1 if peptide has significant MIC activity, 0 otherwise.
    Uses the upper bound of MIC range (>x values map to upper_uM).
    Classification: MIC <= threshold = 1 (active), MIC > threshold or missing = 0 (inactive)
    
    Args:
        lower_uM: Lower bound of concentration in microMolar (str) - not used
        upper_uM: Upper bound of concentration in microMolar (str) - primary classification value
    
    Returns:
        1 if MIC is below threshold (active), 0 otherwise (inactive or missing)
    """
    try:
        # Use upper bound (handles ">x" format which goes to upper_uM)
        if upper_uM:
            mic_value = float(upper_uM)
            return 1 if mic_value <= MIC_THRESHOLD_uM else 0
        # Fall back to lower bound if upper is unavailable
        elif lower_uM:
            mic_value = float(lower_uM)
            return 1 if mic_value <= MIC_THRESHOLD_uM else 0
        # No MIC value available
        else:
            return 0
    except (ValueError, TypeError):
        # If conversion fails, return 0 (no valid activity)
        return 0

def calc_mw(seq: str, nterm: str | None, cterm: str | None) -> float:
    """Calculate peptide molecular weight (Da)."""
    mw = sum(Config.AA_MASS.get(a.upper(), 110.0) for a in (seq or "")) + Config.H2O_MASS
    if nterm and nterm.upper() in Config.NTERM_MASS:
        mw += Config.NTERM_MASS[nterm.upper()]
    if cterm and cterm.upper() in Config.CTERM_MASS:
        mw += Config.CTERM_MASS[cterm.upper()]
    return mw

def get_z_prefix(nterm: str | None) -> str:
    """
    Calculate the number of Z's based on N-terminus.
    Each Z represents a C4 block.
    Examples: C4->Z, C8->ZZ, C12->ZZZ, C16->ZZZZ, C20->ZZZZZ
    
    Args:
        nterm: N-terminus string (e.g., "C12", "C16")
    
    Returns:
        String with appropriate number of Z's
    """
    if not nterm:
        logger.warning("No N-terminus provided, defaulting to C16 (ZZZZ)")
        return "ZZZZ"
    
    match = re.search(r'C(\d+)', nterm.upper())
    if match:
        num = int(match.group(1))
        z_count = num // 4
        if z_count > 0:
            return "Z" * z_count
        else:
            logger.warning(f"Invalid C value {num}, must be multiple of 4. Defaulting to C16 (ZZZZ)")
            return "ZZZZ"
    
    logger.warning(f"Could not parse N-terminus '{nterm}', defaulting to C16 (ZZZZ)")
    return "ZZZZ"

def parse_conc(val: str | None, row_num: int | None = None) -> tuple[str, str]:
    """
    Split concentration string into (lower, upper) strings.
    Handles formats:
    - "value1-value2" -> (value1, value2) [range]
    - "value1->value2" -> (value1, value2) [arrow range]
    - "value±error" -> (value-error, value+error) [mean ± error converted to range]
    - ">value" -> (value, "") [greater than]
    - "<value" -> ("", value) [less than]
    - "<=value" -> ("", value) [less than or equal]
    - "value" -> (value, value) [single value]
    """
    if not val:
        return ("", "")
    # Clean whitespace including tabs
    s = re.sub(r'\s+', ' ', str(val).strip())
    
    # Skip malformed values (like "4.5.5") but allow ranges with decimals
    if s.count('.') > 2 or (s.count('.') > 1 and '±' not in s and '->' not in s and '-' not in s):
        row_info = f" (CSV row {row_num})" if row_num else ""
        logger.warning(f"Malformed concentration value skipped: '{s}'{row_info}")
        return ("", "")
    
    # Handle mean ± error format (e.g., "3.9±1.1" -> lower=2.8, upper=5.0)
    if "±" in s:
        parts = s.split("±", 1)
        try:
            media = float(parts[0].strip())
            error = float(parts[1].strip())
            lower = max(0.0, media - error)  # Don't allow negative concentrations
            upper = media + error
            return (f"{lower:.6g}", f"{upper:.6g}")
        except ValueError:
            row_info = f" (CSV row {row_num})" if row_num else ""
            logger.warning(f"Invalid mean±error format: '{s}'{row_info}")
            return ("", "")
    
    # Handle arrow range format (e.g., "25->50")
    if "->" in s:
        parts = s.split("->", 1)
        return (parts[0].strip(), parts[1].strip())
    
    # Handle less than or equal (e.g., "<=0.25")
    if s.startswith("<="):
        return ("", s[2:].strip())
    
    # Handle greater than (e.g., ">25")
    if s.startswith(">"):
        return ("", s[1:].strip())
    
    # Handle less than (e.g., "<10")
    if s.startswith("<"):
        return ("", s[1:].strip())
    
    # Handle range format (e.g., "6.25-12.5")
    if "-" in s and not s.startswith("-"):
        # Avoid splitting negative numbers
        parts = s.split("-", 1)
        if len(parts) == 2 and parts[0].strip() and parts[1].strip():
            return (parts[0].strip(), parts[1].strip())
    
    # Single value
    return (s, s)

def to_uM(val: str, unit: str | None, mw: float | None, row_num: int | None = None) -> str:
    """Convert one numeric concentration value to µM. Returns '' if not possible."""
    if not val or not unit or not mw:
        return ""
    try:
        num = float(val)
    except ValueError:
        row_info = f" (CSV row {row_num})" if row_num else ""
        logger.warning(f"Invalid concentration value: '{val}'{row_info}")
        return ""
    u = unit.lower()
    if u in ("µm", "um"):
        return f"{num:.6g}"
    if u == "mm":
        return f"{(num * 1000.0):.6g}"
    if u == "nm":
        return f"{(num / 1000.0):.6g}"
    if u == "m":
        return f"{(num * 1_000_000.0):.6g}"
    # mass/volume units require MW (Da ~ g/mol)
    if u in ("µg/ml", "ug/ml", "mg/l"):
        return f"{(num / (mw / 1000.0)):.6g}"
    if u == "g/l":
        return f"{((num * 1000.0) / (mw / 1000.0)):.6g}"
    if u == "ng/ml":
        return f"{((num / 1000.0) / (mw / 1000.0)):.6g}"
    return ""

def split_species(val: str | None) -> tuple[str, str]:
    """
    Split target species into species and strain with a simple rule:
    - skip the first word (usually the genus)
    - if any later token starts with an uppercase letter OR a digit, that token begins the strain
    """
    if not val:
        return ("", "")
    parts = val.split()
    for i, tok in enumerate(parts[1:], start=1):  # skip first word
        t = tok.strip()
        if not t:
            continue
        if t[0].isupper() or t[0].isdigit():
            return (" ".join(parts[:i]), " ".join(parts[i:]))
    return (val, "")

# --- MODIFICADO: cargar mapa {NEW_SEQ -> (curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2)} desde list_min.txt ---
def load_min_map(path: str = None):
    """
    Lee list_min.txt (separado por espacios/tabs) y devuelve un dict:
        key = sequence (por ej., 'ZZZZKLK01')
        value = (curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2) como strings ('' si NA)
    Si hay problemas, devuelve {}.
    """
    if path is None:
        path = Config.MIN_LIST_FILE
        
    mapping = {}
    try:
        with open(path, encoding=Config.CSV_ENCODING) as f:
            header = f.readline().strip().split()
            try:
                idx_seq = header.index("sequence")
                idx_npol = header.index("npol_min")
                idx_curv = header.index("curv_min")
                idx_ph = header.index("pH")
                idx_npol_c0 = header.index("npol_c0")
                idx_npol_c1 = header.index("npol_c1")
                idx_npol_c2 = header.index("npol_c2")
            except ValueError as e:
                logger.warning(f"Missing required columns in {path}: {e}")
                return mapping
                
            for line_num, line in enumerate(f, start=2):
                if not line.strip():
                    continue
                cols = line.strip().split()
                # robustez ante 'NA'
                seq = cols[idx_seq] if len(cols) > idx_seq else ""
                npol = cols[idx_npol] if len(cols) > idx_npol else ""
                curv = cols[idx_curv] if len(cols) > idx_curv else ""
                ph_run = cols[idx_ph] if len(cols) > idx_ph else ""
                npol_c0 = cols[idx_npol_c0] if len(cols) > idx_npol_c0 else ""
                npol_c1 = cols[idx_npol_c1] if len(cols) > idx_npol_c1 else ""
                npol_c2 = cols[idx_npol_c2] if len(cols) > idx_npol_c2 else ""
                
                # Reemplazar 'NA' con cadena vacía
                for val in [npol, curv, ph_run, npol_c0, npol_c1, npol_c2]:
                    if val == "NA":
                        val = ""
                        
                if npol == "NA": npol = ""
                if curv == "NA": curv = ""
                if ph_run == "NA": ph_run = ""
                if npol_c0 == "NA": npol_c0 = ""
                if npol_c1 == "NA": npol_c1 = ""
                if npol_c2 == "NA": npol_c2 = ""
                
                if seq:
                    mapping[seq] = (curv, npol, ph_run, npol_c0, npol_c1, npol_c2)
    except FileNotFoundError:
        logger.warning(f"Min list file not found: {path}")
    except IOError as e:
        logger.error(f"Error reading min list file {path}: {e}")
    return mapping

def run(infile: str | None = None, outfile: str | None = None) -> None:
    """Read activity.csv, compute MW, split/normalize concentrations, split species/strain, write output CSV."""
    if infile is None:
        infile = Config.OUTPUT_ACTIVITY_CSV
    if outfile is None:
        outfile = Config.OUTPUT_NORMALIZED_CSV
        
    logger.info(f"Starting activity normalization: {infile} -> {outfile}")
    
    try:
        # MODIFICADO: cargar join por NEW_SEQ con más columnas
        min_map = load_min_map()
        logger.info(f"Loaded {len(min_map)} entries from min list")

        with open(infile, encoding=Config.CSV_ENCODING) as f:
            r = csv.DictReader(f)
            fieldnames = list(r.fieldnames) + [
                "MW_Da", "NEW_SEQ",
                "lower_concentration", "upper_concentration",
                "lower_uM", "upper_uM",
                "species", "strain",
                # --- COLUMNAS ORIGINALES Y NUEVAS ---
                "curv_min", "npol_min", "ph_run", "npol_c0", "npol_c1", "npol_c2",
                # --- MIC CLASSIFICATION ---
                "has_mic_activity",
            ]
            rows = []
            row_num = 1  # Start at 1, will increment to 2 for first data row
            for row in r:
                row_num += 1  # CSV row number (2 = first data row after header)
                seq = (row.get("SEQUENCE") or "").upper()
                n = row.get("N TERMINUS", "")
                c = row.get("C TERMINUS", "")
                mw = calc_mw(seq, n, c)

                lo, up = parse_conc(row.get("concentration", ""), row_num)
                original_conc = row.get("concentration", "")

                row["MW_Da"] = f"{mw:.2f}"
                z_prefix = get_z_prefix(n)
                new_seq = f"{z_prefix}{seq}{'00' if (c or '').upper() == 'AMD' else '01'}"
                row["NEW_SEQ"] = new_seq
                row["lower_concentration"] = lo
                row["upper_concentration"] = up
                
                # Convert both lower and upper concentrations to µM
                # For ± format: lo=media-error, up=media+error (both are concentrations)
                row["lower_uM"] = to_uM(lo, row.get("unit", ""), mw, row_num)
                row["upper_uM"] = to_uM(up, row.get("unit", ""), mw, row_num)

                sp, st = split_species(row.get("targetSpecies", ""))
                row["species"] = sp
                row["strain"] = st

                # --- Join con curv_min / npol_min / ph_run / npol_c0,c1,c2 por NEW_SEQ ---
                curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2 = ("", "", "", "", "", "")
                if new_seq in min_map:
                    curv_min, npol_min, ph_run, npol_c0, npol_c1, npol_c2 = min_map[new_seq]
                    
                row["curv_min"] = curv_min
                row["npol_min"] = npol_min
                row["ph_run"] = ph_run
                row["npol_c0"] = npol_c0
                row["npol_c1"] = npol_c1
                row["npol_c2"] = npol_c2
                
                # --- Classify MIC activity based on threshold ---
                lower_um = row.get("lower_uM", "")
                upper_um = row.get("upper_uM", "")
                row["has_mic_activity"] = classify_mic_activity(lower_um, upper_um)

                rows.append(row)

            logger.info(f"Processed {len(rows)} rows for normalization")

        with open(outfile, "w", encoding=Config.CSV_ENCODING, newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
            
        logger.info(f"Successfully wrote normalized data to {outfile}")
        
    except FileNotFoundError:
        raise FileProcessingError(f"Input file not found: {infile}", filename=infile)
    except csv.Error as e:
        raise FileProcessingError(f"CSV processing error: {e}", filename=infile)
    except IOError as e:
        raise FileProcessingError(f"File I/O error: {e}", filename=outfile)
    except Exception as e:
        logger.error(f"Unexpected error in normalization: {e}")
        raise

if __name__ == "__main__":
    run()
