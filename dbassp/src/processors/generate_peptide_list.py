import csv
import re
import logging

logger = logging.getLogger("dbaasp_pipeline")

def get_z_count_from_nterm(n_terminus: str) -> int:
    """
    Calculate the number of Z's based on N-terminus.
    Each Z represents a C4 block.
    
    Args:
        n_terminus: N-terminus string (e.g., "C12", "C16")
    
    Returns:
        Number of Z's needed
    """
    match = re.search(r'C(\d+)', n_terminus.upper())
    if match:
        num = int(match.group(1))
        z_count = num // 4
        if z_count > 0:
            return z_count
    
    logger.warning(f"Could not parse N-terminus '{n_terminus}', defaulting to 4 Z's (C16)")
    return 4

def generate_peptide_list(n_terminus):
    """
    Genera archivo peptides_{n_terminus}_list.csv con columnas adicionales
    El número de Z's se calcula dinámicamente basándose en el N-terminus.
    """
    input_file = f"data/input/peptides_{n_terminus}.csv"
    output_file = f"data/output/peptides_{n_terminus}_list.csv"
    
    z_count = get_z_count_from_nterm(n_terminus)
    z_prefix = "Z" * z_count
    
    print(f"Processing {n_terminus}: Using {z_count} Z's ({z_prefix})")
    
    rows = []
    seen_sequences = set()
    
    try:
        with open(input_file, 'r', encoding='utf-8-sig', newline='') as infile:
            reader = csv.DictReader(infile)
            fieldnames = list(reader.fieldnames) + ["total_sequence", "filtered_sequence"]
            
            for row in reader:
                sequence = str(row["SEQUENCE"]).upper()
                c_terminus = str(row["C TERMINUS"]).strip()
                
                if c_terminus.upper() == "AMD":
                    suffix = "00"
                else:
                    suffix = "01"
                    
                total_sequence = f"{z_prefix}{sequence}{suffix}"
                row["total_sequence"] = total_sequence
                
                is_valid = "X" not in total_sequence and total_sequence not in seen_sequences
                
                if is_valid:
                    row["filtered_sequence"] = total_sequence
                    seen_sequences.add(total_sequence)
                else:
                    row["filtered_sequence"] = ""
                
                rows.append(row)
        
        with open(output_file, 'w', encoding='utf-8-sig', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        
        total_sequences = len(rows)
        valid_sequences = sum(1 for row in rows if row["filtered_sequence"])
        sequences_with_x = sum(1 for row in rows if "X" in row["total_sequence"])
        
        print(f"Archivo generado: {output_file}")
        print(f"Total de secuencias originales: {total_sequences}")
        print(f"Secuencias con X (filtradas): {sequences_with_x}")
        print(f"Secuencias únicas válidas: {valid_sequences}")
        print(f"Z prefix usado: {z_prefix} ({z_count} Z's para {n_terminus})")
        
    except FileNotFoundError:
        print(f"Error: No se encontró el archivo {input_file}")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        n_term = sys.argv[1]
        generate_peptide_list(n_term)
    else:
        for n_term in ["C12", "C16"]:
            try:
                generate_peptide_list(n_term)
                print()
            except Exception as e:
                print(f"Error procesando {n_term}: {e}")
