import sqlite3
import requests
import time
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.collectors import lipophilicity as lipo

DB = "dbassp/data/output/dbaasp.sqlite"
conn = sqlite3.connect(DB)
cur = conn.cursor()

def fetch_pubchem_smiles(cid: str) -> str | None:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES,SMILES/JSON"
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=20)
        if resp.ok:
            data = resp.json()
            props = data["PropertyTable"]["Properties"][0]
            return props.get("IsomericSMILES") or props.get("SMILES")
    except Exception:
        pass
    return None

cur.execute("SELECT id, sequence, n_terminus, c_terminus FROM peptides WHERE sequence LIKE '%X%' AND smiles IS NULL")
rows = cur.fetchall()

total = len(rows)
print(f"Found {total} peptides with X and no SMILES.")

updated = 0
for i, (pid, seq, n_term, c_term) in enumerate(rows, 1):
    smiles = None
    try:
        resp = requests.get(f'https://dbaasp.org/peptides/{pid}', headers={'Accept': 'application/json'}, timeout=10)
        if resp.ok:
            data = resp.json()
            
            # Try getting from DBAASP
            s_list = data.get('smiles')
            if s_list and isinstance(s_list, list) and len(s_list) > 0:
                smiles = s_list[0].get('smiles')
            
            # Fallback to PubChem
            if not smiles:
                cid = (data.get("pubChem") or {}).get("cid")
                if cid:
                    smiles = fetch_pubchem_smiles(str(cid))
                    time.sleep(0.3)
                    
            if smiles:
                # Calculate logp / logd
                logp = lipo.calculate_logp(smiles)
                logd = lipo.calculate_logd(
                    smiles, sequence=seq, ph=7.0,
                    nterminus=n_term, cterminus=c_term
                )
                cur.execute("UPDATE peptides SET smiles = ?, logp = ?, logd = ? WHERE id = ?", (smiles, logp, logd, pid))
                updated += 1
                print(f"[{i}/{total}] Updated pid={pid} with SMILES (len {len(smiles)})")
            
        time.sleep(0.5)
    except Exception as e:
        print(f"[{i}/{total}] Error processing pid={pid}: {e}")
        
conn.commit()
conn.close()
print(f"Done! Updated {updated} peptides.")
