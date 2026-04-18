import sqlite3
import glob
import os
import csv
import json
import logging
import time
import requests

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
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            peptide_id INTEGER,
            property_name TEXT,
            property_value TEXT,
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
    for prop in physchem_props:
        name = (prop.get("name") or "").strip()
        val = str(prop.get("value", "")).strip()
        if name and name.upper() != "ID":
            cur.execute("""
                INSERT INTO physchem_properties (peptide_id, property_name, property_value)
                VALUES (?, ?, ?)
            """, (pid, name, val))
            
    # Insert Activity
    target_activities = data.get("targetActivities") or []
    for ta in target_activities:
        ts = (ta.get("targetSpecies") or {}).get("name", "")
        tg = (ta.get("targetGroup") or {}).get("name", "")
        to = (ta.get("targetObject") or {}).get("name", "")
        conc = str(ta.get("concentration", ""))
        measure = (ta.get("activityMeasure") or {}).get("name", "")
        unit = (ta.get("unit") or {}).get("name", "")
        
        cur.execute("""
            INSERT INTO activities (peptide_id, target_species, target_group, target_object, concentration, activity_measure, unit)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (pid, ts, tg, to, conc, measure, unit))
    conn.commit()

def main():
    db_path = "dbaasp_data.sqlite"
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
                        cur.execute('''
                            INSERT INTO peptides (id, complexity, name, n_terminus, sequence, c_terminus, synthesis_type, target_group, target_object)
                            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                        ''', (
                            pid, row.get("COMPLEXITY"), row.get("NAME"), row.get("N TERMINUS"), row.get("SEQUENCE"),
                            row.get("C TERMINUS"), row.get("SYNTHESIS TYPE"), row.get("TARGET GROUP"), row.get("TARGET OBJECT")
                        ))
                    
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
