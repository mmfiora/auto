import sqlite3
import requests

conn = sqlite3.connect('dbassp/data/output/dbaasp.sqlite')
cur = conn.cursor()
cur.execute("SELECT id FROM peptides WHERE sequence LIKE '%X%' AND smiles IS NULL LIMIT 20")
pids = [row[0] for row in cur.fetchall()]

for pid in pids:
    resp = requests.get(f'https://dbaasp.org/peptides/{pid}', headers={'Accept': 'application/json'})
    data = resp.json()
    s = data.get('smiles')
    if s:
        print(f'Peptide {pid} has SMILES in API: {s}')
print('Done checking 20 peptides.')
