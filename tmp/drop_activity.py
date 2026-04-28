import sqlite3

DB = "dbassp/data/output/dbaasp.sqlite"
conn = sqlite3.connect(DB)
cur = conn.cursor()

try:
    cur.execute("ALTER TABLE normalized_activity DROP COLUMN activity")
    conn.commit()
    print("Successfully dropped 'activity' column from 'normalized_activity' table.")
except sqlite3.OperationalError as e:
    print(f"Error (column might already be dropped or SQLite version too old): {e}")

conn.close()
