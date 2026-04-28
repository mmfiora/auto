import sqlite3
import json

DB = "dbassp/data/output/dbaasp.sqlite"

def add_column_back():
    conn = sqlite3.connect(DB)
    cur = conn.cursor()
    
    # 1. Add column back if missing
    cur.execute("PRAGMA table_info(peptides)")
    existing_cols = {r[1] for r in cur.fetchall()}
    if "unusual_residues" not in existing_cols:
        cur.execute("ALTER TABLE peptides ADD COLUMN unusual_residues TEXT")
    
    # 2. Get data from junction table and rebuild JSON
    cur.execute("SELECT peptide_id, position, aa_name FROM peptide_unusual_residues")
    rows = cur.fetchall()
    
    peptide_map = {}
    for pid, pos, aa in rows:
        if pid not in peptide_map:
            peptide_map[pid] = {}
        peptide_map[pid][str(pos)] = aa
        
    # 3. Update peptides table
    update_count = 0
    for pid, unusual_dict in peptide_map.items():
        json_str = json.dumps(unusual_dict)
        cur.execute("UPDATE peptides SET unusual_residues = ? WHERE id = ?", (json_str, pid))
        update_count += 1
        
    conn.commit()
    conn.close()
    print(f"Added 'unusual_residues' column back and populated {update_count} rows.")

if __name__ == '__main__':
    add_column_back()
