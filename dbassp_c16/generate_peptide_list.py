import csv

def generate_peptide_list(n_terminus):
    """
    Genera archivo peptides_{n_terminus}_list.csv con columnas adicionales
    """
    input_file = f"peptides_{n_terminus}.csv"
    output_file = f"peptides_{n_terminus}_list.csv"
    
    # Determinar número de Z según el N terminus
    if n_terminus == "C12":
        z_count = 3  # C12 = 3 grupos de C4
    elif n_terminus == "C16":
        z_count = 4  # C16 = 4 grupos de C4
    else:
        raise ValueError(f"N terminus no soportado: {n_terminus}")
    
    z_prefix = "Z" * z_count
    
    # Leer archivo original y procesar
    rows = []
    seen_sequences = set()
    
    with open(input_file, 'r', encoding='utf-8-sig', newline='') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ["total_sequence", "filtered_sequence"]
        
        for row in reader:
            # Crear total_sequence
            sequence = str(row["SEQUENCE"]).upper()
            c_terminus = str(row["C TERMINUS"]).strip()
            
            # Agregar sufijo basado en C TERMINUS
            if c_terminus == "AMD":
                suffix = "00"
            else:
                suffix = "01"
                
            total_sequence = f"{z_prefix}{sequence}{suffix}"
            row["total_sequence"] = total_sequence
            
            # Determinar si es válida (sin X y no duplicada)
            is_valid = "X" not in total_sequence and total_sequence not in seen_sequences
            
            if is_valid:
                row["filtered_sequence"] = total_sequence
                seen_sequences.add(total_sequence)
            else:
                row["filtered_sequence"] = ""
            
            rows.append(row)
    
    # Escribir archivo de salida
    with open(output_file, 'w', encoding='utf-8-sig', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    
    # Estadísticas
    total_sequences = len(rows)
    valid_sequences = sum(1 for row in rows if row["filtered_sequence"])
    sequences_with_x = sum(1 for row in rows if "X" in row["total_sequence"])
    
    print(f"Archivo generado: {output_file}")
    print(f"Total de secuencias originales: {total_sequences}")
    print(f"Secuencias con X (filtradas): {sequences_with_x}")
    print(f"Secuencias únicas válidas: {valid_sequences}")

if __name__ == "__main__":
    # Generar para C12 y C16
    try:
        generate_peptide_list("C12")
        print()
        generate_peptide_list("C16")
    except Exception as e:
        print(f"Error: {e}")