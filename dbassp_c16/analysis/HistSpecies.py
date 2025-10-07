# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 16:45:01 2025

@author: Maria

Usage:
    python HistSpecies.py                               # Interactive mode
    python HistSpecies.py activity_normalized_C12.csv   # Direct file
    python HistSpecies.py 1                             # By number
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
import sys

# Get all activity_normalized CSV files from data/output directory
output_dir = "data/output"
csv_files = sorted(glob.glob(os.path.join(output_dir, "activity_normalized*.csv")))

if not csv_files:
    print(f"No activity_normalized CSV files found in {output_dir}")
    exit(1)

selected_file = None

# Check if filename was provided as command line argument
if len(sys.argv) > 1:
    arg = sys.argv[1].strip()
    
    # Check if it's a number
    if arg.isdigit():
        idx = int(arg) - 1
        if 0 <= idx < len(csv_files):
            selected_file = csv_files[idx]
            print(f"Seleccionado: {os.path.basename(selected_file)}")
        else:
            print(f"Error: Número fuera de rango. Usa 1-{len(csv_files)}")
            exit(1)
    else:
        # Try to find the file
        potential_file = os.path.join(output_dir, arg)
        if os.path.exists(potential_file):
            selected_file = potential_file
            print(f"Seleccionado: {os.path.basename(selected_file)}")
        elif os.path.exists(potential_file + ".csv"):
            selected_file = potential_file + ".csv"
            print(f"Seleccionado: {os.path.basename(selected_file)}")
        else:
            print(f"Error: Archivo '{arg}' no encontrado en {output_dir}")
            print("\nArchivos disponibles:")
            for i, file in enumerate(csv_files, 1):
                print(f"  {i}. {os.path.basename(file)}")
            exit(1)

# If no argument provided, show interactive menu
if selected_file is None:
    # Display available files
    print("\nArchivos activity_normalized disponibles en data/output/:")
    print("-" * 50)
    for i, file in enumerate(csv_files, 1):
        filename = os.path.basename(file)
        print(f"{i}. {filename}")

    print("-" * 50)

    # Get user selection
    while True:
        try:
            choice = input("\nSelecciona el número del archivo (o escribe el nombre completo): ").strip()
            
            # Check if user entered a number
            if choice.isdigit():
                idx = int(choice) - 1
                if 0 <= idx < len(csv_files):
                    selected_file = csv_files[idx]
                    break
                else:
                    print(f"Error: Selecciona un número entre 1 y {len(csv_files)}")
            else:
                # User entered a filename
                # Try exact match first
                potential_file = os.path.join(output_dir, choice)
                if os.path.exists(potential_file):
                    selected_file = potential_file
                    break
                # Try with .csv extension if not provided
                elif os.path.exists(potential_file + ".csv"):
                    selected_file = potential_file + ".csv"
                    break
                else:
                    print(f"Error: Archivo '{choice}' no encontrado en {output_dir}")
        except KeyboardInterrupt:
            print("\nOperación cancelada")
            exit(0)
        except Exception as e:
            print(f"Error: {e}")

    print(f"\nCargando archivo: {os.path.basename(selected_file)}")

# Load the CSV file
df = pd.read_csv(selected_file)

species_counts=df["species"].value_counts()

plt.figure(figsize=(14, 8))
species_counts.plot(kind="bar")


plt.xlabel("Species", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
plt.xticks(rotation=90)  # etiquetas verticales
plt.tight_layout()
plt.show()