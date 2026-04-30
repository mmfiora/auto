import os
import sys
import sqlite3
import logging

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.collectors import lipophilicity as lipo

logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
log = logging.getLogger(__name__)

DB_PATH = "dbassp/data/output/dbaasp.sqlite"

def main():
    if not os.path.exists(DB_PATH):
        log.error(f"Database not found at {DB_PATH}")
        return

    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    # Fetch all peptides
    cur.execute("SELECT id, sequence, n_terminus, c_terminus, smiles FROM peptides")
    rows = cur.fetchall()
    
    total = len(rows)
    log.info(f"Recalculating logP and logD for {total} peptides...")

    updated_count = 0
    errors_count = 0

    for i, (pid, seq, n_term, c_term, existing_smiles) in enumerate(rows, 1):
        if i % 100 == 0:
            log.info(f"Processed {i}/{total} peptides...")
            
        try:
            smiles = existing_smiles
            
            # For normal peptides (no X), regenerate the SMILES to apply any N-term changes
            if "X" not in seq:
                new_smiles = lipo.sequence_to_smiles(seq, nterminus=n_term, cterminus=c_term)
                if new_smiles:
                    smiles = new_smiles
            
            # If we still don't have a SMILES, we can't calculate logP/logD
            if not smiles:
                continue

            # Recalculate logP and logD
            logp = lipo.calculate_logp(smiles)
            logd = lipo.calculate_logd(smiles, sequence=seq, ph=7.0, nterminus=n_term, cterminus=c_term)

            # Update the database
            cur.execute(
                """UPDATE peptides 
                   SET smiles = ?, logp = ?, logd = ? 
                   WHERE id = ?""",
                (smiles, logp, logd, pid)
            )
            updated_count += 1

        except Exception as e:
            log.warning(f"Failed to recalculate for peptide ID {pid}: {e}")
            errors_count += 1

    conn.commit()
    conn.close()

    log.info(f"Done! Successfully updated {updated_count} peptides.")
    if errors_count > 0:
        log.warning(f"Encountered errors for {errors_count} peptides.")

if __name__ == "__main__":
    main()
