import csv, os, sys
sys.path.insert(0, 'dbassp')
os.chdir('dbassp')
from src.collectors.normalize_activity import load_min_map, get_z_prefix

min_map = load_min_map('data/input/list_min_c16.txt')

with open('data/output/physchem_c16.csv', encoding='utf-8-sig') as f:
    rows = list(csv.DictReader(f))

print(f"{'SEQ':<14}{'C-term':<9}{'our_key(AMD=00)':<22}{'in?':<6}{'alt_key(AMD=01)':<22}{'in?'}")
print('-'*80)
for row in rows:
    seq = (row.get('SEQUENCE') or '').upper()
    n   = row.get('N TERMINUS', '')
    c   = row.get('C TERMINUS', '')
    z   = get_z_prefix(n)
    is_amd = c.upper() == 'AMD'
    our = z + seq + ('00' if is_amd else '01')
    alt = z + seq + ('01' if is_amd else '00')
    our_hit = 'YES' if our in min_map else 'no'
    alt_hit = 'YES' if alt in min_map else 'no'
    if our_hit == 'YES' or alt_hit == 'YES':
        print(f"{seq:<14}{c:<9}{our:<22}{our_hit:<6}{alt:<22}{alt_hit}")
