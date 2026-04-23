import sqlite3

conn = sqlite3.connect("dbassp/data/output/dbaasp.sqlite")
cur = conn.cursor()

print("=== New schema columns ===")
cur.execute("PRAGMA table_info(normalized_activity)")
cols = cur.fetchall()
for r in cols:
    print(r)
new_col_names = {r[1] for r in cols}
assert "c_min_ugml" in new_col_names, "MISSING: c_min_ugml"
assert "c_max_ugml" in new_col_names, "MISSING: c_max_ugml"
assert "conc_gt"    in new_col_names, "MISSING: conc_gt"
assert "activity"   in new_col_names, "MISSING: activity"
print("✓ All 4 new columns present")

print("\n=== Activity distribution ===")
cur.execute("SELECT activity, COUNT(*) FROM normalized_activity GROUP BY activity")
for r in cur.fetchall():
    print(r)

print("\n=== >X rows (conc_gt=1): c_min should be X, c_max NULL ===")
cur.execute("""SELECT concentration, unit, c_min_ugml, c_max_ugml, conc_gt, activity
               FROM normalized_activity WHERE conc_gt=1 LIMIT 10""")
for r in cur.fetchall():
    print(r)

print("\n=== Single value rows: c_min=X, c_max=NULL ===")
cur.execute("""SELECT concentration, unit, c_min_ugml, c_max_ugml, conc_gt, activity
               FROM normalized_activity WHERE conc_gt=0 AND c_max_ugml IS NULL AND c_min_ugml IS NOT NULL LIMIT 10""")
for r in cur.fetchall():
    print(r)

print("\n=== Interval rows: c_min=A, c_max=B ===")
cur.execute("""SELECT concentration, unit, c_min_ugml, c_max_ugml, conc_gt, activity
               FROM normalized_activity WHERE c_max_ugml IS NOT NULL AND conc_gt=0 LIMIT 10""")
for r in cur.fetchall():
    print(r)

print("\n=== Correctness check: active rows must have c_min_ugml <= 32 ===")
cur.execute("SELECT MIN(c_min_ugml), MAX(c_min_ugml) FROM normalized_activity WHERE activity='active'")
mn, mx = cur.fetchone()
print(f"  active c_min range: {mn} .. {mx}")
assert mx is None or mx <= 32, f"FAIL: active row has c_min_ugml={mx} > 32"
print("  ✓ all active rows <= 32 µg/ml")

print("\n=== Correctness check: not-active rows with conc_gt=0 must have c_min_ugml > 32 ===")
cur.execute("SELECT MIN(c_min_ugml), MAX(c_min_ugml) FROM normalized_activity WHERE activity='not active' AND conc_gt=0")
mn2, mx2 = cur.fetchone()
print(f"  not-active (no gt) c_min range: {mn2} .. {mx2}")
assert mn2 is None or mn2 > 32, f"FAIL: not-active row has c_min_ugml={mn2} <= 32"
print("  ✓ all non-gt not-active rows > 32 µg/ml")

print("\n=== Correctness check: conc_gt=1 'unknown' rows must have c_min_ugml <= 32 ===")
cur.execute("SELECT MAX(c_min_ugml) FROM normalized_activity WHERE activity='unknown' AND conc_gt=1")
r = cur.fetchone()
if r[0] is not None:
    assert r[0] <= 32, f"FAIL: unknown gt row has c_min_ugml={r[0]} > 32"
print(f"  max unknown/gt c_min_ugml: {r[0]}  ✓")

conn.close()
print("\nAll checks passed ✓")
