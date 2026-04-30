import os
import glob
import re
import sys
import subprocess

def main():
    # Find project root (2 levels up from scripts/collection/fetch_all_to_csv.py)
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    
    # Find all peptide input files
    peptide_files = glob.glob(os.path.join(project_root, "data", "input", "peptides_*.csv"))
    if not peptide_files:
        print(f"Error: No peptides files found in {os.path.join(project_root, 'data', 'input')}")
        return 1
    
    pattern = re.compile(r'peptides_([A-Za-z]*\d+)\.csv$', re.IGNORECASE)
    
    ntermini = []
    for f in sorted(peptide_files):
        basename = os.path.basename(f)
        match = pattern.match(basename)
        if match:
            nt = match.group(1)
            ntermini.append(nt)
            
    if not ntermini:
        print("Error: Could not extract any N-termini from file names.")
        return 1
        
    print(f"Found N-termini to process: {', '.join(ntermini)}")
    
    failed = []
    main_api_path = os.path.join(project_root, "main_api.py")
    for nt in ntermini:
        print(f"\n{'='*50}\nProcessing N-terminus: {nt}\n{'='*50}")
        # Run main_api.py via subprocess to ensure clean state and avoid module import conflicts
        result = subprocess.run([sys.executable, main_api_path, "--nterminus", nt], cwd=project_root)
        if result.returncode != 0:
            print(f"Error processing {nt}")
            failed.append(nt)
            
    if failed:
        print(f"\nCompleted with errors for: {', '.join(failed)}")
        return 1
        
    print("\nSuccessfully fetched all data to CSVs for all N-termini!")
    return 0

if __name__ == "__main__":
    sys.exit(main())
