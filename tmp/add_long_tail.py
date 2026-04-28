import sqlite3
import re

DB = "dbassp/data/output/dbaasp.sqlite"

conn = sqlite3.connect(DB)
cur = conn.cursor()

# Add the column (ignore if it already exists)
try:
    cur.execute("ALTER TABLE peptides ADD COLUMN long_tail INTEGER")
    print("Added column: long_tail")
except Exception as e:
    print(f"Column may already exist: {e}")

# Populate: extract the number after 'C' in n_terminus
cur.execute("SELECT id, n_terminus FROM peptides WHERE n_terminus IS NOT NULL")
rows = cur.fetchall()

updated = 0
for pid, n_term in rows:
    m = re.fullmatch(r'C(\d+)', n_term.strip())
    if m:
        cur.execute("UPDATE peptides SET long_tail = ? WHERE id = ?", (int(m.group(1)), pid))
        updated += 1

conn.commit()
print(f"Updated {updated} rows")

# Verify
cur.execute("SELECT n_terminus, long_tail FROM peptides WHERE long_tail IS NOT NULL LIMIT 10")
print("\nSample:")
for r in cur.fetchall():
    print(f"  n_terminus={r[0]!r:6}  long_tail={r[1]}")

cur.execute("SELECT COUNT(*) FROM peptides WHERE n_terminus IS NOT NULL AND long_tail IS NULL")
missed = cur.fetchone()[0]
print(f"\nRows with n_terminus but no long_tail (unmatched): {missed}")

conn.close()
