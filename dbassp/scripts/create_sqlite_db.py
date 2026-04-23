import sqlite3
import glob
import os
import sys
import csv
import logging
import time
import requests

# Add the parent directory to sys.path to import from src
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from src.collectors import lipophilicity
except ImportError:
    lipophilicity = None
    logging.warning("Failed to import lipophilicity module. SMILES/logP/logD will not be calculated.")

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')

def fetch_peptide(pid):
    url = f"https://dbaasp.org/peptides/{pid}"
    headers = {"Accept": "application/json"}
    try:
        resp = requests.get(url, headers=headers, timeout=20)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        logging.error(f"Error fetching {pid}: {e}")
        return None

def setup_db(conn):
    cur = conn.cursor()
    cur.execute('''
        CREATE TABLE IF NOT EXISTS peptides (
            id              INTEGER PRIMARY KEY,
            complexity      TEXT,
            name            TEXT,
            n_terminus      TEXT,
            sequence        TEXT,
            c_terminus      TEXT,
            synthesis_type  TEXT,
            target_group    TEXT,
            target_object   TEXT,
            smiles          TEXT,
            logp            REAL,
            logd            REAL
        )
    ''')
    cur.execute('''
        CREATE TABLE IF NOT EXISTS physchem_properties_aa (
            peptide_id                              INTEGER PRIMARY KEY,
            normalized_hydrophobic_moment           REAL,
            normalized_hydrophobicity               REAL,
            net_charge                              REAL,
            isoelectric_point                       REAL,
            penetration_depth                       REAL,
            tilt_angle                              REAL,
            disordered_conformation_propensity      REAL,
            linear_moment                           REAL,
            propensity_to_in_vitro_aggregation      REAL,
            angle_subtended_by_the_hydrophobic_residues REAL,
            amphiphilicity_index                    REAL,
            propensity_to_ppii_coil                 REAL,
            FOREIGN KEY(peptide_id) REFERENCES peptides(id)
        )
    ''')
    cur.execute('''
        CREATE TABLE IF NOT EXISTS normalized_activity (
            id               INTEGER PRIMARY KEY AUTOINCREMENT,
            peptide_id       INTEGER NOT NULL REFERENCES peptides(id),
            target_species   TEXT,
            species          TEXT,
            strain           TEXT,
            target_group     TEXT,
            target_object    TEXT,
            concentration    TEXT,
            c_min_uM         REAL,
            c_max_uM         REAL,
            c_min_ugml       REAL,
            c_max_ugml       REAL,
            conc_gt          INTEGER,
            activity         TEXT,
            activity_measure TEXT,
            unit             TEXT,
            ph               TEXT,
            ionic_strength   TEXT,
            salt_type        TEXT,
            medium           TEXT,
            cfu              TEXT,
            note             TEXT,
            reference        TEXT
        )
    ''')
    cur.execute('CREATE INDEX IF NOT EXISTS idx_physchem_peptide ON physchem_properties_aa(peptide_id)')
    cur.execute('CREATE INDEX IF NOT EXISTS idx_normact_peptide  ON normalized_activity(peptide_id)')
    conn.commit()

# ── Concentration / unit helpers (mirrors migrate_db.py) ──────────────────────
import re

_PHYSCHEM_MAP = {
    'Normalized Hydrophobic Moment': 'normalized_hydrophobic_moment',
    'Normalized Hydrophobicity':     'normalized_hydrophobicity',
    'Net Charge':                    'net_charge',
    'Isoelectric Point':             'isoelectric_point',
    'Penetration Depth':             'penetration_depth',
    'Tilt Angle':                    'tilt_angle',
    'Disordered Conformation Propensity': 'disordered_conformation_propensity',
    'Linear Moment':                 'linear_moment',
    'Propensity to in vitro Aggregation': 'propensity_to_in_vitro_aggregation',
    'Angle Subtended by the Hydrophobic Residues': 'angle_subtended_by_the_hydrophobic_residues',
    'Amphiphilicity Index':          'amphiphilicity_index',
    'Propensity to PPII coil':       'propensity_to_ppii_coil',
}

def _parse_conc(val):
    """Return (lower_raw, upper_raw, conc_gt). conc_gt=1 if form was '>X' or '>=X'."""
    if not val:
        return ("", "", 0)
    s = re.sub(r'\s+', ' ', str(val).strip())
    # mean ± error
    if '\u00b1' in s:
        parts = s.split('\u00b1', 1)
        try:
            mu  = float(parts[0].strip())
            err = float(parts[1].strip())
            return (f"{max(0.0, mu-err):.6g}", f"{mu+err:.6g}", 0)
        except ValueError:
            return ("", "", 0)
    if '->' in s:
        a, b = s.split('->', 1); return (a.strip(), b.strip(), 0)
    # >X / >=X  →  c_min = X, c_max = NULL, conc_gt = 1
    if s.startswith('>='):  return (s[2:].strip(), "", 1)
    if s.startswith('>'):   return (s[1:].strip(), "", 1)
    # <X / <=X  →  c_min = NULL, c_max = X
    if s.startswith('<='): return ("", s[2:].strip(), 0)
    if s.startswith('<'):  return ("", s[1:].strip(), 0)
    if '-' in s and not s.startswith('-'):
        parts = s.split('-', 1)
        if parts[0].strip() and parts[1].strip():
            return (parts[0].strip(), parts[1].strip(), 0)
    # Single value  →  c_min = X, c_max = NULL
    return (s, "", 0)

def _to_uM(val, unit, mw):
    if not val: return None
    try: num = float(val)
    except ValueError: return None
    u = (unit or '').strip().lower().replace('\ufffd', 'µ').replace('\xb5', 'µ').replace('?', 'µ')
    if u in ('µm', 'um', '\u03bcm'): return round(num, 6)
    if u == 'mm': return round(num * 1e3, 6)
    if u == 'nm': return round(num / 1e3, 6)
    if u == 'm':  return round(num * 1e6, 6)
    if u in ('µg/ml', 'ug/ml', '\u03bcg/ml', 'mg/l', 'mg/ml'):
        return round(num / (mw / 1e3), 6) if mw else None
    if u in ('g/l', 'g/ml'):
        return round((num * 1e3) / (mw / 1e3), 6) if mw else None
    if u in ('ng/ml', 'ng/l'):
        return round((num / 1e3) / (mw / 1e3), 6) if mw else None
    return None

def _to_ugml(val, unit, mw):
    if not val: return None
    try: num = float(val)
    except ValueError: return None
    u = (unit or '').strip().lower().replace('\ufffd', 'µ').replace('\xb5', 'µ').replace('?', 'µ')
    if u in ('µg/ml', 'ug/ml', '\u03bcg/ml', 'mg/l'): return round(num, 6)
    if u == 'mg/ml': return round(num * 1e3, 6)
    if u in ('g/l', 'g/ml'): return round(num * 1e3, 6)
    if u in ('ng/ml', 'ng/l'): return round(num / 1e3, 6)
    if not mw: return None
    if u in ('µm', 'um', '\u03bcm'): return round(num * (mw / 1e3), 6)
    if u == 'mm': return round(num * 1e3 * (mw / 1e3), 6)
    if u == 'nm': return round(num / 1e3 * (mw / 1e3), 6)
    if u == 'm':  return round(num * 1e6 * (mw / 1e3), 6)
    return None

def _classify_activity(c_min_ugml, conc_gt):
    """Return 'active', 'not active', or 'unknown' (threshold: 32 µg/ml)."""
    if c_min_ugml is None:
        return 'unknown'
    if conc_gt:
        return 'not active' if c_min_ugml > 32 else 'unknown'
    return 'active' if c_min_ugml <= 32 else 'not active'

def _split_species(val):
    if not val: return ("", "")
    parts = val.split()
    for i, tok in enumerate(parts[1:], start=1):
        if tok.strip() and (tok[0].isupper() or tok[0].isdigit()):
            return (" ".join(parts[:i]), " ".join(parts[i:]))
    return (val, "")

# ─────────────────────────────────────────────────────────────────────────────

def process_api_data(conn, pid, data, mw=None):
    """Insert physchem (AA-only) and normalized activity rows for peptide pid."""
    if not data:
        return
    cur = conn.cursor()

    # ── physchem_properties_aa ────────────────────────────────────────────────
    phys_data = {"peptide_id": pid}
    for prop in data.get("physicoChemicalProperties") or []:
        name = (prop.get("name") or "").strip()
        val  = str(prop.get("value", "")).strip()
        if name in _PHYSCHEM_MAP:
            phys_data[_PHYSCHEM_MAP[name]] = val

    if len(phys_data) > 1:
        cols   = list(phys_data.keys())
        values = list(phys_data.values())
        cur.execute(
            f"INSERT INTO physchem_properties_aa ({', '.join(cols)}) "
            f"VALUES ({', '.join(['?']*len(cols))})",
            values
        )

    # ── normalized_activity ───────────────────────────────────────────────────
    for ta in data.get("targetActivities") or []:
        ts      = (ta.get("targetSpecies") or {}).get("name", "")
        tg      = (ta.get("targetGroup")   or {}).get("name", "")
        tobj    = (ta.get("targetObject")  or {}).get("name", "")
        conc    = str(ta.get("concentration", ""))
        measure = (ta.get("activityMeasure") or {}).get("name", "")
        unit    = (ta.get("unit")   or {}).get("name", "")
        ph      = str(ta.get("ph")           or "")
        ionic   = str(ta.get("ionicStrength") or "")
        salt    = str(ta.get("saltType")      or "")
        medium  = str((ta.get("medium") or {}).get("name", "") if isinstance(ta.get("medium"), dict) else ta.get("medium") or "")
        cfu     = str(ta.get("cfu")  or "")
        note    = str(ta.get("note") or "")
        ref     = str(ta.get("reference") or "")

        lo, hi, gt = _parse_conc(conc)
        c_min     = _to_uM(lo, unit, mw)
        c_max     = _to_uM(hi, unit, mw)
        c_min_ugl = _to_ugml(lo, unit, mw)
        c_max_ugl = _to_ugml(hi, unit, mw)
        activity  = _classify_activity(c_min_ugl, gt)
        species, strain = _split_species(ts)

        cur.execute("""
            INSERT INTO normalized_activity
               (peptide_id, target_species, species, strain,
                target_group, target_object, concentration,
                c_min_uM, c_max_uM, c_min_ugml, c_max_ugml, conc_gt, activity,
                activity_measure, unit,
                ph, ionic_strength, salt_type, medium, cfu, note, reference)
            VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """, (pid, ts, species, strain, tg, tobj, conc,
               c_min, c_max, c_min_ugl, c_max_ugl, gt, activity,
               measure, unit, ph, ionic, salt, medium, cfu, note, ref))

    conn.commit()

def main():
    db_path = "data/output/dbaasp.sqlite"
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    setup_db(conn)

    csv_files = glob.glob("data/input/peptides_*.csv")
    if not csv_files:
        logging.error("No peptides_*.csv found in data/input/")
        return

    cur = conn.cursor()
    for file_path in csv_files:
        logging.info(f"Processing {file_path}...")
        try:
            with open(file_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    id_raw = row.get("ID", "")
                    if not id_raw:
                        continue
                    try:
                        pid = int(id_raw.split("_")[-1])
                    except ValueError:
                        continue
                    
                    # Insert peptide base info (with inline smiles/logp/logd)
                    cur.execute("SELECT 1 FROM peptides WHERE id = ?", (pid,))
                    if not cur.fetchone():
                        n_term = row.get("N TERMINUS")
                        seq    = row.get("SEQUENCE")
                        c_term = row.get("C TERMINUS")

                        smiles = logp = logd = None
                        if lipophilicity is not None and seq:
                            smiles = lipophilicity.sequence_to_smiles(seq, nterminus=n_term, cterminus=c_term)
                            if smiles:
                                logp = lipophilicity.calculate_logp(smiles)
                                logd = lipophilicity.calculate_logd(smiles, sequence=seq, ph=7.0, nterminus=n_term, cterminus=c_term)

                        cur.execute('''
                            INSERT INTO peptides
                               (id, complexity, name, n_terminus, sequence, c_terminus,
                                synthesis_type, target_group, target_object,
                                smiles, logp, logd)
                            VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
                        ''', (
                            pid, row.get("COMPLEXITY"), row.get("NAME"), n_term, seq,
                            c_term, row.get("SYNTHESIS TYPE"),
                            row.get("TARGET GROUP"), row.get("TARGET OBJECT"),
                            smiles, logp, logd
                        ))
                    
                    # Check if physchem data already exists (resume-safe)
                    cur.execute("SELECT 1 FROM physchem_properties_aa WHERE peptide_id = ? LIMIT 1", (pid,))
                    if cur.fetchone():
                        logging.debug(f"API data for peptide {pid} already processed. Skipping.")
                        continue

                    # Compute approximate MW for concentration normalisation
                    seq_   = row.get("SEQUENCE") or ""
                    n_term_ = row.get("N TERMINUS") or ""
                    AA_MASS = {'A':71.08,'R':156.19,'N':114.10,'D':115.09,'C':103.15,
                               'E':129.12,'Q':128.13,'G':57.05,'H':137.14,'I':113.16,
                               'L':113.16,'K':128.17,'M':131.20,'F':147.18,'P':97.12,
                               'S':87.08,'T':101.11,'W':186.21,'Y':163.18,'V':99.13,
                               'Z':56.10,'X':110.0}
                    NTERM_MASS = {'C16':239.2,'C14':211.2,'C12':183.2,'C10':155.2,
                                  'C8':127.2,'C18':267.2,'C20':295.2,'C6':99.2,'C4':71.2}
                    mw = sum(AA_MASS.get(a.upper(), 110.0) for a in seq_) + 18.02
                    mw += NTERM_MASS.get(n_term_.upper(), 0.0)

                    logging.info(f"Fetching API data for peptide {pid}...")
                    data = fetch_peptide(pid)
                    if data:
                        process_api_data(conn, pid, data, mw=mw)
                    time.sleep(0.5)  # Rate limiting
        except Exception as e:
            logging.error(f"Failed to process {file_path}: {e}")

    conn.close()
    logging.info(f"Database creation complete successfully! Saved to {db_path}")

if __name__ == "__main__":
    main()
