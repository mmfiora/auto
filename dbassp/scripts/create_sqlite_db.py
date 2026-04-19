import sqlite3
import glob
import os
import sys
import csv
import json
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
            id INTEGER PRIMARY KEY,
            complexity TEXT,
            name TEXT,
            n_terminus TEXT,
            sequence TEXT,
            c_terminus TEXT,
            synthesis_type TEXT,
            target_group TEXT,
            target_object TEXT
        )
    ''')
    cur.execute('''
        CREATE TABLE IF NOT EXISTS physchem_properties (
            peptide_id INTEGER PRIMARY KEY,
            normalized_hydrophobic_moment REAL,
            normalized_hydrophobicity REAL,
            net_charge REAL,
            isoelectric_point REAL,
            penetration_depth REAL,
            tilt_angle REAL,
            disordered_conformation_propensity REAL,
            linear_moment REAL,
            propensity_to_in_vitro_aggregation REAL,
            angle_subtended_by_the_hydrophobic_residues REAL,
            amphiphilicity_index REAL,
            propensity_to_ppii_coil REAL,
            FOREIGN KEY(peptide_id) REFERENCES peptides(id)
        )
    ''')
    cur.execute('''
        CREATE TABLE IF NOT EXISTS activities (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            peptide_id INTEGER,
            target_species TEXT,
            target_group TEXT,
            target_object TEXT,
            concentration TEXT,
            activity_measure TEXT,
            unit TEXT,
            ph TEXT,
            ionic_strength TEXT,
            salt_type TEXT,
            medium TEXT,
            cfu TEXT,
            note TEXT,
            reference TEXT,
            FOREIGN KEY(peptide_id) REFERENCES peptides(id)
        )
    ''')
    cur.execute('''
        CREATE TABLE IF NOT EXISTS lipophilicity (
            peptide_id INTEGER PRIMARY KEY,
            smiles TEXT,
            logp REAL,
            logd REAL,
            FOREIGN KEY(peptide_id) REFERENCES peptides(id)
        )
    ''')
    cur.execute('CREATE INDEX IF NOT EXISTS idx_physchem_peptide ON physchem_properties(peptide_id)')
    cur.execute('CREATE INDEX IF NOT EXISTS idx_activities_peptide ON activities(peptide_id)')
    conn.commit()

def process_api_data(conn, pid, data):
    if not data:
        return
    cur = conn.cursor()
    # Insert Physchem
    physchem_props = data.get("physicoChemicalProperties") or []
    PHYSCHEM_MAP = {
        'Normalized Hydrophobic Moment': 'normalized_hydrophobic_moment',
        'Normalized Hydrophobicity': 'normalized_hydrophobicity',
        'Net Charge': 'net_charge',
        'Isoelectric Point': 'isoelectric_point',
        'Penetration Depth': 'penetration_depth',
        'Tilt Angle': 'tilt_angle',
        'Disordered Conformation Propensity': 'disordered_conformation_propensity',
        'Linear Moment': 'linear_moment',
        'Propensity to in vitro Aggregation': 'propensity_to_in_vitro_aggregation',
        'Angle Subtended by the Hydrophobic Residues': 'angle_subtended_by_the_hydrophobic_residues',
        'Amphiphilicity Index': 'amphiphilicity_index',
        'Propensity to PPII coil': 'propensity_to_ppii_coil'
    }
    
    phys_data = {"peptide_id": pid}
    for prop in physchem_props:
        name = (prop.get("name") or "").strip()
        val = str(prop.get("value", "")).strip()
        if name in PHYSCHEM_MAP:
            phys_data[PHYSCHEM_MAP[name]] = val
    
    if len(phys_data) > 1:
        cols = ", ".join(phys_data.keys())
        placeholders = ", ".join(["?"] * len(phys_data))
        cur.execute(f"INSERT INTO physchem_properties ({cols}) VALUES ({placeholders})", list(phys_data.values()))
            
    # Insert Activity
    target_activities = data.get("targetActivities") or []
    for ta in target_activities:
        ts = (ta.get("targetSpecies") or {}).get("name", "")
        tg = (ta.get("targetGroup") or {}).get("name", "")
        to = (ta.get("targetObject") or {}).get("name", "")
        conc = str(ta.get("concentration", ""))
        measure = (ta.get("activityMeasure") or {}).get("name", "")
        unit = (ta.get("unit") or {}).get("name", "")
        ph = str((ta.get("ph") or "") if ta.get("ph") is not None else "")
        ionic = str((ta.get("ionicStrength") or "") if ta.get("ionicStrength") is not None else "")
        salt = str((ta.get("saltType") or "") if ta.get("saltType") is not None else "")
        medium = str((ta.get("medium") or "") if ta.get("medium") is not None else "")
        cfu = str((ta.get("cfu") or "") if ta.get("cfu") is not None else "")
        note = str((ta.get("note") or "") if ta.get("note") is not None else "")
        ref = str((ta.get("reference") or "") if ta.get("reference") is not None else "")
        
        cur.execute("""
            INSERT INTO activities (peptide_id, target_species, target_group, target_object, concentration, activity_measure, unit, ph, ionic_strength, salt_type, medium, cfu, note, reference)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (pid, ts, tg, to, conc, measure, unit, ph, ionic, salt, medium, cfu, note, ref))
    conn.commit()

def main():
    db_path = "data/output/dbaasp_data.sqlite"
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
                    
                    # Insert peptide base info
                    cur.execute("SELECT 1 FROM peptides WHERE id = ?", (pid,))
                    if not cur.fetchone():
                        n_term = row.get("N TERMINUS")
                        seq = row.get("SEQUENCE")
                        c_term = row.get("C TERMINUS")
                        
                        cur.execute('''
                            INSERT INTO peptides (id, complexity, name, n_terminus, sequence, c_terminus, synthesis_type, target_group, target_object)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (
                            pid, row.get("COMPLEXITY"), row.get("NAME"), n_term, seq,
                            c_term, row.get("SYNTHESIS TYPE"), row.get("TARGET GROUP"), row.get("TARGET OBJECT")
                        ))
                        
                        # Calculate and insert lipophilicity
                        if lipophilicity is not None and seq:
                            smiles = lipophilicity.sequence_to_smiles(seq, nterminus=n_term, cterminus=c_term)
                            if smiles:
                                logp = lipophilicity.calculate_logp(smiles)
                                logd = lipophilicity.calculate_logd(smiles, sequence=seq, ph=7.0, nterminus=n_term, cterminus=c_term)
                                cur.execute('''
                                    INSERT INTO lipophilicity (peptide_id, smiles, logp, logd)
                                    VALUES (?, ?, ?, ?)
                                ''', (pid, smiles, logp, logd))
                    
                    # Check if physchem data already exists (so we don't refetch on resume)
                    cur.execute("SELECT 1 FROM physchem_properties WHERE peptide_id = ? LIMIT 1", (pid,))
                    if cur.fetchone():
                        logging.debug(f"API data for peptide {pid} already processed. Skipping.")
                        continue
                        
                    logging.info(f"Fetching API data for peptide {pid}...")
                    data = fetch_peptide(pid)
                    if data:
                        process_api_data(conn, pid, data)
                    time.sleep(0.5)  # Rate limiting
        except Exception as e:
            logging.error(f"Failed to process {file_path}: {e}")

    conn.close()
    logging.info(f"Database creation complete successfully! Saved to {db_path}")

if __name__ == "__main__":
    main()
