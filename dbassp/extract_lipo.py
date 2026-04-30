import sqlite3
import csv
import os

db_path = "data/output/dbaasp.sqlite"
conn = sqlite3.connect(db_path)
cur = conn.cursor()

cur.execute("SELECT DISTINCT n_terminus FROM peptides")
ntermini = [r[0] for r in cur.fetchall() if r[0]]

for nt in ntermini:
    cur.execute("SELECT id, smiles, logp, logd FROM peptides WHERE n_terminus = ?", (nt,))
    rows = cur.fetchall()
    
    out_file = f"data/output/lipophilicity_{nt}.csv"
    with open(out_file, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(["Peptide ID", "SMILES", "logP", "logD"])
        for r in rows:
            writer.writerow(r)
    print(f"Exported {out_file} with {len(rows)} rows.")

conn.close()
