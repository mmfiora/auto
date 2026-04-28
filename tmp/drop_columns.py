import sqlite3

DB = "dbassp/data/output/dbaasp.sqlite"

# Columns still to drop (c_min_uM was already dropped successfully)
cols_to_drop = ["c_max_uM", "curv_min", "npol_min", "ph_run", "npol_c0", "npol_c1", "npol_c2"]

conn = sqlite3.connect(DB)
cur = conn.cursor()

# 1. Drop the view that references these columns
cur.execute("DROP VIEW IF EXISTS peptides_activity_view")
print("Dropped view: peptides_activity_view")

# 2. Drop the columns
for col in cols_to_drop:
    cur.execute(f"ALTER TABLE normalized_activity DROP COLUMN {col}")
    print(f"Dropped column: {col}")

# 3. Recreate the view without the dropped columns
cur.execute("""
CREATE VIEW peptides_activity_view AS
SELECT ROW_NUMBER() OVER () AS row_id, p.*, na.species, na.strain
FROM peptides p
JOIN normalized_activity na ON p.id = na.peptide_id
""")
print("Recreated view: peptides_activity_view")

conn.commit()

# Verify
print("\n--- normalized_activity columns ---")
cur.execute("PRAGMA table_info(normalized_activity)")
for r in cur.fetchall():
    print(f"  {r[1]}")

print("\n--- peptides_activity_view columns ---")
cur.execute("PRAGMA table_info(peptides_activity_view)")
for r in cur.fetchall():
    print(f"  {r[1]}")

conn.close()
