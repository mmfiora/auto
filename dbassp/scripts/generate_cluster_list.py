import csv
import os
import sys
import argparse

def generate_cluster_list(input_file: str, output_file: str):
    if not os.path.exists(input_file):
        print(f"Error: Could not find {input_file}")
        return False
        
    cluster_strings = set()
    total_processed = 0
    filtered_out_x = 0
    
    with open(input_file, 'r', encoding='utf-8-sig', newline='') as infile:
        reader = csv.DictReader(infile)
        if 'SEQUENCE' not in reader.fieldnames or 'N TERMINUS' not in reader.fieldnames or 'C TERMINUS' not in reader.fieldnames:
            print(f"Error: Missing required columns (SEQUENCE, N TERMINUS, C TERMINUS) in {input_file}.")
            return False
            
        for row in reader:
            total_processed += 1
            sequence = str(row['SEQUENCE']).strip().upper()
            n_term = str(row['N TERMINUS']).strip().upper()
            c_term = str(row['C TERMINUS']).strip().upper()
            
            # Filter unusual amino acids (X)
            if 'X' in sequence:
                filtered_out_x += 1
                continue
                
            # Determine prefix
            if n_term == 'C16':
                prefix = 'ZZZZ'
            elif n_term == 'C12':
                prefix = 'ZZZ'
            else:
                # If it's not C16 or C12, omit from this list
                continue
                
            # Determine suffix
            if c_term == 'AMD':
                suffix = '00'
            else:
                suffix = '01'
                
            # Combine to get final string (e.g. ZZZZKAAK00)
            final_string = f"{prefix}{sequence}{suffix}"
            cluster_strings.add(final_string)
            
    # Write to text file as a raw list
    with open(output_file, 'w', encoding='utf-8') as outfile:
        for s in sorted(cluster_strings):
            outfile.write(f"{s}\n")
            
    print(f"Processed: {input_file}")
    print(f"  Total sequences parsed: {total_processed}")
    print(f"  Filtered out (contained 'X'): {filtered_out_x}")
    print(f"  Unique valid lists saved: {len(cluster_strings)}")
    print(f"  Saved to: {output_file}\n")
    return True

def main():
    parser = argparse.ArgumentParser(description="Create a plain-text list of peptide sequences for the cluster.")
    parser.add_argument('input_files', nargs='*', help="Path to input CSV database extractions.")
    
    args = parser.parse_args()
    
    # If no files are specified, default to the standard ones
    if not args.input_files:
        files_to_process = [
            "data/input/peptides_C12.csv",
            "data/input/peptides_C16.csv"
        ]
        print(f"No inputs provided. Defaulting to: {', '.join(files_to_process)}\n")
    else:
        files_to_process = args.input_files
        
    success = False
    for filepath in files_to_process:
        if os.path.exists(filepath):
            filename = os.path.basename(filepath)
            name, ext = os.path.splitext(filename)
            # Create a simple .txt for the cluster
            output_filepath = os.path.join(os.path.dirname(filepath), f"{name}_cluster_list.txt")
            if generate_cluster_list(filepath, output_filepath):
                success = True
        else:
            print(f"Warning: {filepath} not found.")
            
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
