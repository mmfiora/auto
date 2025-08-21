# main.py
import time
import physchem
import activity
import normalize_activity

def main():
    t0 = time.perf_counter()

    print("Running physchem ...")
    physchem.run()

    print("Running activity ...")
    activity.run()

    print("Normalizing activity ...")
    normalize_activity.run(infile="activity.csv", outfile="activity_normalized.csv")

    t1 = time.perf_counter()
    print(f"Done. Total elapsed: {t1 - t0:.2f} s")

if __name__ == "__main__":
    main()

