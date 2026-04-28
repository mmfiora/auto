import sqlite3

DB = "dbassp/data/output/dbaasp.sqlite"
conn = sqlite3.connect(DB)
cur = conn.cursor()
print(f'SQLite version: {sqlite3.sqlite_version}')

print('\n--- unusual_amino_acids ---')
cur.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='unusual_amino_acids'")
res = cur.fetchone()
if res: print(res[0])

print('\n--- peptides ---')
cur.execute("SELECT sql FROM sqlite_master WHERE type='table' AND name='peptides'")
res = cur.fetchone()
if res: print(res[0])

conn.close()
