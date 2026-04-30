import sqlite3
import glob
import os
import sys
import csv
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')

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
            logd            REAL,
            molecular_weight REAL,
            total_charge     REAL,
            long_tail        TEXT,
            unusual_aminoacid TEXT -- JSON format {"1": "DAB", ...}
        )
    ''')
    cur.execute('''
        CREATE TABLE IF NOT EXISTS peptide_unusual_aminoacid (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            peptide_id      INTEGER,
            aa_name         TEXT,
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
            activity_measure TEXT,
            unit             TEXT,
            ph               TEXT,
            ionic_strength   TEXT,
            salt_type        TEXT,
            medium           TEXT,
            cfu              TEXT,
            note             TEXT,
            reference        TEXT,
            c_min_ugml       REAL,
            c_max_ugml       REAL,
            conc_gt          INTEGER
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
        CREATE TABLE IF NOT EXISTS unusual_aminoacids (
            id              INTEGER PRIMARY KEY AUTOINCREMENT,
            name            TEXT UNIQUE,
            description     TEXT
        )
    ''')
    cur.execute('CREATE INDEX IF NOT EXISTS idx_physchem_peptide ON physchem_properties_aa(peptide_id)')
    cur.execute('CREATE INDEX IF NOT EXISTS idx_normact_peptide  ON normalized_activity(peptide_id)')
    cur.execute('CREATE INDEX IF NOT EXISTS idx_unusual_peptide  ON peptide_unusual_aminoacid(peptide_id)')
    conn.commit()

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

def to_float_or_null(val):
    if not val: return None
    try: return float(val)
    except ValueError: return None

def process_peptides(conn):
    cur = conn.cursor()
    csv_files = glob.glob("data/input/peptides_*.csv")
    for file_path in csv_files:
        logging.info(f"Populating peptides from {file_path}...")
        try:
            with open(file_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    id_raw = row.get("ID", "")
                    if not id_raw: continue
                    try: pid = int(id_raw.split("_")[-1])
                    except ValueError: continue
                    
                    cur.execute("SELECT 1 FROM peptides WHERE id = ?", (pid,))
                    if cur.fetchone(): continue

                    cur.execute('''
                        INSERT INTO peptides
                           (id, complexity, name, n_terminus, sequence, c_terminus,
                            synthesis_type, target_group, target_object)
                        VALUES (?,?,?,?,?,?,?,?,?)
                    ''', (
                        pid, row.get("COMPLEXITY"), row.get("NAME"), row.get("N TERMINUS"), row.get("SEQUENCE"),
                        row.get("C TERMINUS"), row.get("SYNTHESIS TYPE"),
                        row.get("TARGET GROUP"), row.get("TARGET OBJECT")
                    ))
        except Exception as e:
            logging.error(f"Failed to process {file_path}: {e}")
    conn.commit()

def process_physchem(conn):
    cur = conn.cursor()
    csv_files = glob.glob("data/output/intrinsic_properties_*.csv")
    for file_path in csv_files:
        logging.info(f"Populating intrinsic properties from {file_path}...")
        try:
            with open(file_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    pid_raw = row.get("Peptide ID", "")
                    if not pid_raw: continue
                    try: pid = int(pid_raw)
                    except ValueError: continue
                    
                    # Update peptides table with intrinsic metrics
                    cur.execute('''
                        UPDATE peptides SET 
                            smiles=?, logp=?, logd=?, 
                            molecular_weight=?, total_charge=?, long_tail=?
                        WHERE id=?
                    ''', (
                        row.get("SMILES"), to_float_or_null(row.get("logP")), to_float_or_null(row.get("logD")),
                        to_float_or_null(row.get("molecular_weight")), 
                        to_float_or_null(row.get("total_charge")), 
                        row.get("long_tail"),
                        pid
                    ))
                    # Note: I put SMILES as logP accidentally in the above line's target, fixing it.
                    # Correcting: cur.execute(..., (row.get("SMILES"), ...))
                    
                    # Insert into physchem_properties_aa
                    cur.execute("SELECT 1 FROM physchem_properties_aa WHERE peptide_id = ?", (pid,))
                    if cur.fetchone(): continue
                    
                    phys_data = {"peptide_id": pid}
                    for csv_col, db_col in _PHYSCHEM_MAP.items():
                        if csv_col in row and row[csv_col]:
                            phys_data[db_col] = row[csv_col]
                    
                    if len(phys_data) > 1:
                        cols = list(phys_data.keys())
                        values = list(phys_data.values())
                        cur.execute(
                            f"INSERT INTO physchem_properties_aa ({', '.join(cols)}) VALUES ({', '.join(['?']*len(cols))})",
                            values
                        )
        except Exception as e:
            logging.error(f"Failed to process {file_path}: {e}")
    conn.commit()

import json

def process_activity(conn):
    cur = conn.cursor()
    csv_files = glob.glob("data/output/activity_normalized_*.csv")
    for file_path in csv_files:
        logging.info(f"Populating normalized activity from {file_path}...")
        try:
            with open(file_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    pid_raw = row.get("Peptide ID", "")
                    if not pid_raw: continue
                    try: pid = int(pid_raw)
                    except ValueError: continue
                    
                    unusual_map_raw = row.get("Unusual Amino Acids Map", "{}")
                    if unusual_map_raw and unusual_map_raw != "{}":
                        cur.execute("UPDATE peptides SET unusual_aminoacid=? WHERE id=?", (unusual_map_raw, pid))
                        
                        # Populate peptide_unusual_aminoacid table
                        try:
                            unusual_map = json.loads(unusual_map_raw)
                            for pos, aa_name in unusual_map.items():
                                # Avoid duplicates for the same peptide
                                cur.execute("SELECT 1 FROM peptide_unusual_aminoacid WHERE peptide_id=? AND aa_name=?", (pid, aa_name))
                                if not cur.fetchone():
                                    cur.execute("INSERT INTO peptide_unusual_aminoacid (peptide_id, aa_name) VALUES (?,?)", (pid, aa_name))
                        except json.JSONDecodeError:
                            pass
                        
                    cur.execute("""
                        INSERT INTO normalized_activity
                           (peptide_id, target_species, species, strain,
                            target_group, target_object, concentration,
                            activity_measure, unit,
                            ph, ionic_strength, salt_type, medium, cfu, note, reference,
                            c_min_ugml, c_max_ugml, conc_gt)
                        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                    """, (
                        pid, row.get("targetSpecies", ""), row.get("species", ""), row.get("strain", ""),
                        row.get("targetGroup", ""), row.get("targetObject", ""), row.get("concentration", ""),
                        row.get("activityMeasureGroup", ""), row.get("unit", ""),
                        row.get("ph", ""), row.get("ionicStrength", ""), row.get("saltType", ""),
                        row.get("medium", ""), row.get("cfu", ""), row.get("note", ""), row.get("reference", ""),
                        to_float_or_null(row.get("lower_ugml")), to_float_or_null(row.get("upper_ugml")),
                        int(row.get("conc_gt", 0)) if row.get("conc_gt") else 0
                    ))
        except Exception as e:
            logging.error(f"Failed to process {file_path}: {e}")
    conn.commit()

def main():
    db_path = "data/output/dbaasp.sqlite"
    
    # Optional: Delete existing database to ensure a clean build from CSVs
    if os.path.exists(db_path):
        logging.info(f"Removing existing database {db_path}...")
        os.remove(db_path)
        
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    setup_db(conn)

    process_peptides(conn)
    process_physchem(conn)
    process_activity(conn)

    conn.close()
    logging.info(f"Database creation complete! Saved to {db_path}")

if __name__ == "__main__":
    main()

