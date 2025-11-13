import csv
import logging
import math
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from src.core import common
from src.core.config import Config
from src.core.exceptions import APIError, FileProcessingError

logger = logging.getLogger("dbaasp_pipeline")

N_TERMINUS_MODIFICATIONS = {
    "C12": "CCCCCCCCCCCC(=O)",
    "C16": "CCCCCCCCCCCCCCCC(=O)",
}

# Amino acid pKa values for ionizable side chains
# Acidic amino acids: side chain COOH groups
ACIDIC_AA_PKA = {
    'D': 3.9,   # Aspartate
    'C': 8.37,  # Cysteine
    'E': 4.07,  # Glutamate
    'Y': 10.5,  # Tyrosine
}

# Basic amino acids: side chain NH2/NH3+ groups
# pKb values converted to pKa using: pKa = 14 - pKb
BASIC_AA_PKA = {
    'R': 12.48,  # Arginine (pKb=1.52)
    'H': 6.04,   # Histidine (pKb=7.96)
    'K': 10.54,  # Lysine (pKb=3.46)
    'O': 10.50,  # Ornithine (pKb=3.5)
    'B': 10.27,  # Diaminobutyric acid (pKb=3.73)
    'J': 9.43,   # Diaminopropionic acid (pKb=4.57)
}

# Terminal group pKa values
N_TERMINUS_PKA = 9.0   # NH3+ group
C_TERMINUS_PKA = 3.5   # COOH group

def get_modification_smiles(nterminus: str) -> str | None:
    """Get SMILES for N-terminal modification (C12, C16, etc.)."""
    return N_TERMINUS_MODIFICATIONS.get(nterminus)

def sequence_to_smiles(sequence: str, nterminus: str | None = None, cterminus: str | None = None) -> str | None:
    """
    Convert peptide sequence to SMILES string.
    Optionally append N-terminal modification (C12, C16).
    Handle C-terminal modification (AMD = amidated).
    
    Args:
        sequence: Peptide sequence
        nterminus: N-terminal modification (e.g., "C12", "C16")
        cterminus: C-terminal modification (e.g., "AMD" for amidated, None/empty for free COOH)
    
    Returns:
        SMILES string with appropriate terminal modifications
    """
    try:
        mol = AllChem.MolFromSequence(sequence)
        if mol is None:
            logger.warning(f"Failed to generate SMILES for sequence: {sequence}")
            return None
        
        # Handle C-terminus amidation (AMD)
        # RDKit's MolFromSequence creates COOH by default, we need to convert to CONH2
        if cterminus and cterminus.upper() == "AMD":
            # Convert C-terminus COOH to CONH2 (amide)
            # This is done by editing the molecule to replace OH with NH2
            mol = Chem.RWMol(mol)
            
            # Find the C-terminal carboxyl oxygen (OH)
            # Look for C(=O)O pattern at the end
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1:  # OH group
                    # Check if connected to carbonyl carbon
                    for bond in atom.GetBonds():
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetSymbol() == 'C':
                            # Check if this carbon has a double-bonded oxygen (C=O)
                            has_carbonyl = False
                            for nbond in neighbor.GetBonds():
                                other = nbond.GetOtherAtom(neighbor)
                                if other.GetSymbol() == 'O' and nbond.GetBondType() == Chem.BondType.DOUBLE:
                                    has_carbonyl = True
                                    break
                            
                            if has_carbonyl:
                                # Replace OH with NH2
                                atom_idx = atom.GetIdx()
                                mol.ReplaceAtom(atom_idx, Chem.Atom(7))  # 7 = Nitrogen
                                mol.GetAtomWithIdx(atom_idx).SetNumExplicitHs(2)  # NH2
                                break
                    break
            
            mol = mol.GetMol()
        
        smiles = Chem.MolToSmiles(mol)
        
        # Handle N-terminus modification
        if nterminus:
            mod_smiles = get_modification_smiles(nterminus)
            if mod_smiles:
                smiles = f"{mod_smiles}[N+]({smiles})[H]"
        
        return smiles
    except Exception as e:
        logger.warning(f"Error generating SMILES for sequence {sequence}: {e}")
        return None

def calculate_logp(smiles: str) -> float | None:
    """Calculate logP (partition coefficient) from SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None
        
        logp = Descriptors.MolLogP(mol)
        return logp
    except Exception as e:
        logger.warning(f"Error calculating logP for SMILES {smiles}: {e}")
        return None

def count_ionizable_groups(sequence: str) -> dict:
    """
    Count ionizable groups in peptide sequence.
    Returns dict with counts of acidic and basic residues.
    """
    acidic_count = {}
    basic_count = {}
    
    for aa in sequence.upper():
        if aa in ACIDIC_AA_PKA:
            acidic_count[aa] = acidic_count.get(aa, 0) + 1
        elif aa in BASIC_AA_PKA:
            basic_count[aa] = basic_count.get(aa, 0) + 1
    
    return {
        'acidic': acidic_count,
        'basic': basic_count
    }

def calculate_logd(smiles: str, sequence: str = None, ph: float = 7.4, 
                   nterminus: str = None, cterminus: str = None) -> float | None:
    """
    Calculate logD (distribution coefficient) at specified pH.
    logD is logP adjusted for ionization at given pH using Henderson-Hasselbalch equation.
    
    Args:
        smiles: SMILES string of the molecule
        sequence: Peptide sequence (needed for ionization correction)
        ph: pH value (default 7.4 for physiological pH)
        nterminus: N-terminal modification (e.g., "C12", "C16", None for free NH2)
        cterminus: C-terminal modification (e.g., "AMD" for amidated, None for free COOH)
    
    Returns:
        logD value corrected for ionization, or None if calculation fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            return None
        
        # Start with logP (partition coefficient of neutral form)
        logp = Descriptors.MolLogP(mol)
        
        # If no sequence provided, return logP (no ionization correction)
        if not sequence:
            logger.debug("No sequence provided for logD calculation, returning logP")
            return logp
        
        # Count ionizable groups
        ionizable = count_ionizable_groups(sequence)
        
        # Calculate ionization correction using Henderson-Hasselbalch
        correction = 0.0
        
        # Acidic groups: COOH -> COO- (lose proton at high pH)
        # At pH > pKa, group is deprotonated (charged)
        # Correction: -log10(1 + 10^(pH - pKa))
        for aa, count in ionizable['acidic'].items():
            pka = ACIDIC_AA_PKA[aa]
            ionization_factor = 1 + 10**(ph - pka)
            correction -= count * math.log10(ionization_factor)
        
        # Basic groups: NH3+ -> NH2 (lose proton at high pH)
        # At pH < pKa, group is protonated (charged)
        # Correction: -log10(1 + 10^(pKa - pH))
        for aa, count in ionizable['basic'].items():
            pka = BASIC_AA_PKA[aa]
            ionization_factor = 1 + 10**(pka - ph)
            correction -= count * math.log10(ionization_factor)
        
        # N-terminus correction (if not modified)
        # Modified N-terminus (C12, C16, etc.) is not ionizable
        has_free_n_terminus = not (nterminus and nterminus in N_TERMINUS_MODIFICATIONS)
        if has_free_n_terminus:
            # N-terminus: NH3+ (basic)
            ionization_factor = 1 + 10**(N_TERMINUS_PKA - ph)
            correction -= math.log10(ionization_factor)
        
        # C-terminus correction (if not amidated)
        # AMD (amidated) C-terminus is CONH2 (neutral, not ionizable)
        # Only free COOH C-terminus is ionizable
        has_free_c_terminus = not (cterminus and cterminus.upper() == "AMD")
        if has_free_c_terminus:
            # C-terminus: COOH (acidic)
            ionization_factor = 1 + 10**(ph - C_TERMINUS_PKA)
            correction -= math.log10(ionization_factor)
        
        logd = logp + correction
        
        logger.debug(f"logP={logp:.3f}, correction={correction:.3f}, logD={logd:.3f} at pH {ph}")
        
        return logd
        
    except Exception as e:
        logger.warning(f"Error calculating logD for SMILES {smiles}: {e}")
        return None

def run():
    logger.info("Starting lipophilicity data collection")
    
    try:
        nterminus = Config.get_nterminus()
        logger.info(f"Processing lipophilicity for N-terminus: {nterminus}")
        
        ids = common.load_ids()
        if not ids:
            logger.warning(f"No IDs found in {Config.INPUT_PEPTIDES_CSV}")
            return
        
        logger.info(f"Processing {len(ids)} peptides for lipophilicity data")
        
        data = []
        failed_peptides = []
        
        for i, pid in enumerate(ids, start=1):
            print(f"[{i}/{len(ids)}] Fetching lipophilicity for peptide {pid} ...", end="", flush=True)
            try:
                peptide_data = common.fetch(pid)
                
                sequence = peptide_data.get("sequence", "").strip()
                if not sequence:
                    logger.warning(f"No sequence found for peptide {pid}")
                    failed_peptides.append(pid)
                    print(" fail (no sequence)")
                    continue
                
                # Extract C-terminus information
                cterminus_data = peptide_data.get("cTerminus") or {}
                cterminus = (cterminus_data.get("name") or "").strip()
                
                # Generate SMILES with proper terminal modifications
                smiles = sequence_to_smiles(sequence, nterminus=nterminus, cterminus=cterminus)
                if not smiles:
                    logger.warning(f"Failed to generate SMILES for peptide {pid}")
                    failed_peptides.append(pid)
                    print(" fail (SMILES generation)")
                    continue
                
                logp = calculate_logp(smiles)
                logd = calculate_logd(smiles, sequence=sequence, ph=7.0, 
                                     nterminus=nterminus, cterminus=cterminus)
                
                data.append({
                    "Peptide ID": str(peptide_data.get("id", "")),
                    "N TERMINUS": nterminus,
                    "SEQUENCE": sequence,
                    "SMILES": smiles,
                    "logP": logp if logp is not None else "",
                    "logD": logd if logd is not None else "",
                })
                print(" ok")
                
            except APIError as e:
                logger.error(f"API error for peptide {pid}: {e}")
                failed_peptides.append(pid)
                print(f" fail (API error)")
            except Exception as e:
                logger.error(f"Unexpected error for peptide {pid}: {e}")
                failed_peptides.append(pid)
                print(f" fail (unexpected error)")
        
        if failed_peptides:
            logger.warning(f"Failed to process {len(failed_peptides)} peptides: {failed_peptides}")
        
        if not data:
            logger.error("No data calculated; nothing to write.")
            return
        
        logger.info(f"Successfully calculated lipophilicity for {len(data)} peptides")
        
        header = ["Peptide ID", "N TERMINUS", "SEQUENCE", "SMILES", "logP", "logD"]
        
        try:
            with open(Config.OUTPUT_LIPOPHILICITY_CSV, "w", newline="", encoding=Config.CSV_ENCODING) as f:
                w = csv.DictWriter(f, fieldnames=header)
                w.writeheader()
                w.writerows(data)
            
            logger.info(f"Successfully wrote lipophilicity data to {Config.OUTPUT_LIPOPHILICITY_CSV}")
            
        except IOError as e:
            raise FileProcessingError(f"Error writing lipophilicity file: {e}", filename=Config.OUTPUT_LIPOPHILICITY_CSV)
    
    except (FileProcessingError, APIError) as e:
        logger.error(f"Lipophilicity collection failed: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in lipophilicity collection: {e}")
        raise
