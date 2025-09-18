# activity.py
import csv
import common
from config import Config

def get_activities(peptide_json):
    a = peptide_json.get("targetActivities")
    if a is None:
        a = peptide_json.get("activityAgainstTargetSpecies")
    return a if isinstance(a, list) else []

def collect_activity_keys(all_peptides):
    exclude = {"id", "activity", "activityMeasureValue"}
    order, seen = [], set()
    for d in all_peptides:
        for a in get_activities(d):
            for k in a.keys():
                if not k or k.lower() in exclude:
                    continue
                if k not in seen:
                    seen.add(k)
                    order.append(k)
    return order

def flatten_value(v):
    if v is None:
        return ""
    if isinstance(v, dict):
        return v.get("name", str(v))
    return v

def run():
    ids = common.load_ids()
    if not ids:
        print(f"No IDs found in {Config.INPUT_PEPTIDES_CSV}")
        return

    peptides = []
    for i, pid in enumerate(ids, start=1):
        print(f"[{i}/{len(ids)}] Fetching activity for peptide {pid} ...", end="", flush=True)
        try:
            peptides.append(common.fetch(pid))
            print(" ok")
        except Exception as e:
            print(f" fail ({e})")

    if not peptides:
        print("No data fetched; nothing to write.")
        return

    activity_cols = collect_activity_keys(peptides)
    header = ["Peptide ID", "N TERMINUS", "SEQUENCE", "C TERMINUS"] + activity_cols

    with open(Config.OUTPUT_ACTIVITY_CSV, "w", newline="", encoding=Config.CSV_ENCODING) as f:
        w = csv.DictWriter(f, fieldnames=header)
        w.writeheader()

        for d in peptides:
            base = {
                "Peptide ID": str(d.get("id", "")),
                "N TERMINUS": (d.get("nTerminus") or {}).get("name", ""),
                "SEQUENCE": d.get("sequence", ""),
                "C TERMINUS": (d.get("cTerminus") or {}).get("name", ""),
            }
            acts = get_activities(d)
            if not acts:
                w.writerow(base)
                continue
            for a in acts:
                row = dict(base)
                for k in activity_cols:
                    row[k] = flatten_value(a.get(k))
                w.writerow(row)


