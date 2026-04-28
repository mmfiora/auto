import sqlite3
import requests

conn = sqlite3.connect('dbassp/data/output/dbaasp.sqlite')
cur = conn.cursor()
cur.execute("SELECT id FROM peptides WHERE sequence LIKE '%X%' AND smiles IS NULL LIMIT 1")
pid = cur.fetchone()[0]

resp = requests.get(f'https://dbaasp.org/peptides/{pid}', headers={'Accept': 'application/json'})
data = resp.json()

print('Keys:', list(data.keys()))
if 'smiles' in data:
    print('SMILES:', data['smiles'])
else:
    print('No direct SMILES field')

print('PubChem data:', data.get('pubChem'))
