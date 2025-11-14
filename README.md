# DBAASP Peptide Analysis Pipeline

A Python pipeline for analyzing antimicrobial peptides (AMPs) from the DBAASP (Database of Antimicrobial Activity and Structure of Peptides) database. This tool collects and processes physicochemical properties and biological activity data for peptides with different N-terminal modifications (C12, C16, etc.).

## Overview

Automates collection and analysis of antimicrobial peptides (AMPs) from the DBAASP database:
- Fetches peptide data via DBAASP API
- Extracts physicochemical properties, biological activity, and lipophilicity metrics
- Normalizes activity data for comparative analysis across species
- Supports multiple N-terminal modifications (C12, C16, etc.)

## Project Structure

```
dbassp_c16/
├── data/
│   ├── input/           # Input CSV files with peptide sequences
│   └── output/          # Generated analysis output files
├── src/
│   ├── collectors/      # Data collection modules
│   │   ├── activity.py          # Biological activity data collection
│   │   ├── normalize_activity.py # Activity normalization
│   │   └── physchem.py          # Physicochemical properties collection
│   ├── core/           # Core utilities and configuration
│   │   ├── common.py           # Shared functions and API interaction
│   │   ├── config.py           # Configuration management
│   │   └── exceptions.py       # Custom exceptions
│   └── processors/     # Data processing modules
│       ├── generate_peptide_list.py # Peptide list generation
│       └── unified_results.py       # Results consolidation
├── analysis/           # Analysis scripts
│   ├── CorrelationMatrix.py
│   └── HistSpecies.py
└── main_api.py         # Main entry point
```

## Requirements

- Python 3.12+
- Dependencies: See `requirements.txt`
- Internet connection (for DBAASP API access)

## Installation

1. Clone or download the project
2. Create and activate a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. **Prepare Input**: Place peptide data in `data/input/peptides_{NTERMINUS}.csv` with columns:
   - `Peptide ID` or `ID`
   - `N TERMINUS` (e.g., "C16")
   - `SEQUENCE`

2. **Run Pipeline**:
   ```bash
   python main_api.py
   ```
   Auto-detects N-terminus from filename, fetches data from DBAASP API, extracts properties, and generates output files

### Input Files

- **list_min_{NTERMINUS}.txt** - Text file with peptide sequences and MOLT (Molecular Theory)-calculated minima. Includes free energy, aggregation number (npol), and preferred curvature (curv_min) obtained from cluster runs. Values are merged into unified_results. 
- **peptides_{NTERMINUS}.csv** - Peptide dataset with ID, sequence, and N-terminal modifications

### Output Files

- **physchem_{NTERMINUS}.csv** - Physicochemical properties
- **activity_{NTERMINUS}.csv** - Biological activity data
- **lipophilicity_{NTERMINUS}.csv** - Lipophilicity metrics (logP, logD, free energy, aggregation number, curvature)
- **activity_normalized_{NTERMINUS}.csv** - Normalized activity metrics
- **unified_results_{NTERMINUS}.csv** - Combined analysis results

## Configuration

Configuration is managed via `src/core/config.py`. Key settings can be overridden with environment variables:

```bash
export DBAASP_API_URL="https://dbaasp.org/peptides/{id}"
export DBAASP_TIMEOUT="20"
export LOG_LEVEL="DEBUG"  # Set to DEBUG for verbose output
```

File paths are automatically set based on the N-terminus modification (C12, C16, etc.)

## Extracted Data

**Physicochemical** (DBAASP): Molecular weight, pI, charge, hydrophobicity, aromaticity, instability index, aliphatic index, amino acid composition

**Biological Activity**: Target species, activity measurements, MIC, antibacterial/antifungal potency

**Self-Assembly** (MOLT): Free energy, aggregation number (npol), self-assembly curvature

## Troubleshooting

**Missing Input Files**: Ensure peptide CSV files exist in `data/input/` with pattern `peptides_{NTERMINUS}.csv`

**API Connection Issues**:
- Check internet connection and DBAASP API accessibility
- Increase `DBAASP_TIMEOUT` if requests time out

**Data Validation Errors**:
- Verify CSV files contain required columns: `Peptide ID` or `ID`, `N TERMINUS`, `SEQUENCE`
- Check UTF-8 encoding and valid peptide ID format

**Debug Mode**:
```bash
export LOG_LEVEL="DEBUG"
python main_api.py
```

## Support

For issues: check `logs/pipeline.log`, enable DEBUG logging (`export LOG_LEVEL="DEBUG"`), or review source code documentation
