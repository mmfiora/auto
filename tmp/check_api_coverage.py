import sqlite3
import requests

conn = sqlite3.connect('dbassp/data/output/dbaasp.sqlite')
cur = conn.cursor()
cur.execute("SELECT id FROM peptides WHERE sequence LIKE '%X%'")
pids = [row[0] for row in cur.fetchall()]

has_smiles_count = 0
for pid in pids:
    resp = requests.get(f'https://dbaasp.org/peptides/{pid}', headers={'Accept': 'application/json'})
    if not resp.ok: continue
    data = resp.json()
    s = data.get('smiles')
    if s and len(s) > 0 and s[0].get('smiles'):
        has_smiles_count += 1
        
print(f'Out of {len(pids)} peptides with X, {has_smiles_count} have SMILES in DBAASP API.')
