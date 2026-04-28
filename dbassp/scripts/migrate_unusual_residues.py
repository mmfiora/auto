import sqlite3
import json

DB = "dbassp/data/output/dbaasp.sqlite"

def migrate():
    conn = sqlite3.connect(DB)
    cur = conn.cursor()
    
    # 1. Create junction table
    cur.execute("""
        CREATE TABLE IF NOT EXISTS peptide_unusual_residues (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            peptide_id INTEGER NOT NULL,
            position INTEGER NOT NULL,
            aa_name TEXT NOT NULL,
            FOREIGN KEY (peptide_id) REFERENCES peptides (id),
            FOREIGN KEY (aa_name) REFERENCES unusual_amino_acids (name)
        )
    """)
    
    # 2. Extract data from JSON and insert into junction table
    cur.execute("SELECT id, unusual_residues FROM peptides WHERE unusual_residues IS NOT NULL")
    rows = cur.fetchall()
    
    insert_count = 0
    for peptide_id, unusual_residues_json in rows:
        try:
            data = json.loads(unusual_residues_json)
            for pos_str, aa_name in data.items():
                position = int(pos_str)
                cur.execute(
                    "INSERT INTO peptide_unusual_residues (peptide_id, position, aa_name) VALUES (?, ?, ?)",
                    (peptide_id, position, aa_name)
                )
                insert_count += 1
        except json.JSONDecodeError:
            print(f"Warning: Could not parse JSON for peptide {peptide_id}: {unusual_residues_json}")
    
    # 3. Drop old column
    try:
        cur.execute("ALTER TABLE peptides DROP COLUMN unusual_residues")
        print("Successfully dropped 'unusual_residues' column from 'peptides' table.")
    except sqlite3.OperationalError as e:
        print(f"Notice: Could not drop column (maybe already dropped or SQLite version too old): {e}")

    conn.commit()
    conn.close()
    
    print(f"Migration completed. Inserted {insert_count} junction records.")

if __name__ == '__main__':
    migrate()
