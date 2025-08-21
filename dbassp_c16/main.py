# main.py
import os
import time
import csv
from get_peptide_card import get_peptide_card_info
from extract_activity import write_activity_csv
from extract_physchem import write_physchem_csv

def main():
    ids = [
        51, 838, 860, 862, 863, 864, 1124, 1125, 6197, 6198, 6199,
        7322, 7342, 7348, 7542, 7580, 7927, 7973, 8348, 8473, 8506,
        9257, 9259, 9261, 9263, 9277, 9278, 9279, 9280, 9281, 9282,
        9286, 9288, 9291, 9292, 9294, 9295, 9296, 9297, 9298, 9299,
        9300, 9301, 9302, 9303, 9304, 9305, 9306, 9307, 9309, 9745,
        9746, 10361, 10362, 10675, 10678, 10680, 11017, 11018, 11033,
        11034, 11035, 11150, 11151, 11155, 11156, 11157, 11179, 11182,
        12103, 12106, 12348, 12362, 14155, 14172, 14284, 14288, 15268,
        16365, 16853, 16888, 16890, 17440, 17441, 17442, 17443, 17444,
        17445, 17446, 17447, 17736, 17742, 17784, 17786, 17846, 18097,
        18902, 22614, 22617, 22620, 22623, 22626, 22629, 23366, 23904,
        23907
    ]

    txt_files = []
    per_id_times = []  # lista de dicts: {"Peptide ID": id, "t_s": segundos}

    t_download_start = time.perf_counter()
    for pid in ids:
        print(f"Procesando ID {pid}...")
        t0 = time.perf_counter()
        get_peptide_card_info(pid)
        t1 = time.perf_counter()

        txt_path = f"peptide_{pid}.txt"
        if os.path.isfile(txt_path):
            txt_files.append(txt_path)
            per_id_times.append({
                "Peptide ID": pid,
                "t_s": round(t1 - t0, 3)
            })
        else:
            per_id_times.append({
                "Peptide ID": pid,
                "t_s": None
            })
    t_download_end = time.perf_counter()

    if txt_files:
        t_act_start = time.perf_counter()
        n_act = write_activity_csv(txt_files, "activity.csv")
        t_act_end = time.perf_counter()

        t_pc_start = time.perf_counter()
        n_pc = write_physchem_csv(txt_files, "physchem.csv")
        t_pc_end = time.perf_counter()

        # Borrar los txt después de procesarlos
        for txt in txt_files:
            try:
                os.remove(txt)
                print(f"Archivo {txt} borrado.")
            except Exception as e:
                print(f"No se pudo borrar {txt}: {e}")

        # Guardar timings por ID a CSV
        with open("timings.csv", "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=["Peptide ID", "t_s"])
            w.writeheader()
            for row in per_id_times:
                w.writerow(row)

        # Resumen en consola
        descargados = [r for r in per_id_times if r["t_s"] is not None]
        n_ok = len(descargados)
        t_descarga_total = t_download_end - t_download_start
        t_act = t_act_end - t_act_start
        t_pc = t_pc_end - t_pc_start
        t_total = t_descarga_total + t_act + t_pc
        prom = (sum(r["t_s"] for r in descargados) / n_ok) if n_ok else 0.0

        print("==== Resumen de tiempos ====")
        print(f"IDs totales: {len(ids)}")
        print(f"IDs descargados OK: {n_ok}")
        print(f"Tiempo descarga total: {t_descarga_total:.2f} s")
        print(f"Tiempo escribir activity.csv: {t_act:.2f} s (filas: {n_act})")
        print(f"Tiempo escribir physchem.csv: {t_pc:.2f} s (filas: {n_pc})")
        print(f"Promedio por péptido (descarga): {prom:.2f} s")
        print(f"Tiempo total: {t_total:.2f} s")
    else:
        print("No hay archivos peptide_<id>.txt para exportar.")

if __name__ == "__main__":
    main()

