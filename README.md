# DBAASP Peptide Analysis Pipeline

A comprehensive Python pipeline for analyzing antimicrobial peptides (AMPs) from the DBAASP (Database of Antimicrobial Activity and Structure of Peptides) database. This tool collects, processes, and analyzes physicochemical properties and biological activity data for peptides with different N-terminal modifications (C12, C16, etc.).

## Overview

This pipeline automates the collection and analysis of peptide data through several key stages:

1. **Data Collection** - Fetches peptide data from the DBAASP API
2. **Physicochemical Analysis** - Extracts molecular properties (mass, charge, hydrophobicity, etc.)
3. **Activity Analysis** - Extracts and processes biological activity information
4. **Lipophilicity Analysis** - Calculates logP and logD from peptide sequences
5. **Normalization** - Normalizes activity values for comparative analysis
6. **Unified Results** - Combines all data into comprehensive analysis files

## Features

- **Multi-terminal Support**: Handles peptides with different N-terminal modifications (C12, C16, etc.)
- **API Integration**: Direct integration with DBAASP API for real-time data fetching
- **Comprehensive Analysis**: Extracts physicochemical properties per peptide
- **Lipophilicity Calculation**: Computes logP and logD for peptides with N-terminal modifications
- **Activity Tracking**: Captures detailed activity metrics across multiple target species
- **Data Validation**: Built-in validation and error handling for data integrity
- **Configurable**: Environment-based configuration with sensible defaults
- **Logging**: Detailed logging for debugging and monitoring pipeline execution

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
│   │   ├── physchem.py          # Physicochemical properties extraction
│   │   └── lipophilicity.py      # Lipophilicity (logP/logD) calculation
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

### Basic Workflow

1. **Prepare Input Data**: Create a CSV file with peptide sequences in `data/input/` with format:
   ```
   peptides_{NTERMINUS}.csv
   ```
   Example: `peptides_C16.csv`

   Required columns:
   - `Peptide ID` or `ID` - Unique peptide identifier
   - `N TERMINUS` - N-terminal modification (e.g., "C16")
   - `SEQUENCE` - Amino acid sequence
   - Other optional peptide properties

2. **Run the Pipeline**:
   ```bash
   python main_api.py
   ```

    The pipeline will:
    - Auto-detect the N-terminal modification from the input filename
    - Fetch peptide data from DBAASP API
    - Extract physicochemical properties
    - Calculate lipophilicity (logP/logD) for peptides
    - Extract activity information
    - Generate normalized and unified output files

### Output Files

The pipeline generates several output CSV files:

- **physchem_{NTERMINUS}.csv** - Physicochemical properties for each peptide
- **activity_{NTERMINUS}.csv** - Biological activity data
- **activity_normalized_{NTERMINUS}.csv** - Normalized activity metrics
- **lipophilicity_{NTERMINUS}.csv** - Lipophilicity descriptors (logP, logD) and SMILES
- **unified_results_{NTERMINUS}.csv** - Combined analysis of all properties and activities

## Configuration

Configuration is managed via `dbassp_c16/src/core/config.py` with environment variable overrides:

### Key Configuration Options

```python
# API Settings
DBAASP_API_URL = "https://dbaasp.org/peptides/{id}"
DBAASP_TIMEOUT = 20  # seconds

# File Paths (set automatically based on N-terminus)
INPUT_PEPTIDES_CSV = "data/input/peptides_{NTERMINUS}.csv"
OUTPUT_PHYSCHEM_CSV = "data/output/physchem_{NTERMINUS}.csv"
OUTPUT_ACTIVITY_CSV = "data/output/activity_{NTERMINUS}.csv"
OUTPUT_LIPOPHILICITY_CSV = "data/output/lipophilicity_{NTERMINUS}.csv"
OUTPUT_NORMALIZED_CSV = "data/output/activity_normalized_{NTERMINUS}.csv"
OUTPUT_UNIFIED_CSV = "data/output/unified_results_{NTERMINUS}.csv"

# Logging
LOG_LEVEL = "INFO"  # Set to DEBUG for verbose output
```

### Environment Variables

Override defaults by setting environment variables:

```bash
export DBAASP_API_URL="https://custom.dbaasp.org/peptides/{id}"
export DBAASP_TIMEOUT="30"
export LOG_LEVEL="DEBUG"
python main_api.py
```

## Physicochemical Properties Calculated

The pipeline calculates the following properties for each peptide:

- Molecular weight
- Isoelectric point (pI)
- Charge at different pH values
- Hydrophobicity indices
- Aromaticity
- Instability index
- Aliphatic index
- Amino acid composition
- And more...

## Biological Activity Metrics

The pipeline captures:

- Target organisms/species
- Activity measurements and units
- Minimum inhibitory concentration (MIC)
- Antibacterial/antifungal potency
- Cross-species activity comparisons

## Troubleshooting

### Missing Input Files
Ensure you have peptide CSV files in `data/input/` with the pattern `peptides_{NTERMINUS}.csv`

### API Connection Errors
- Check internet connection
- Verify DBAASP API is accessible at the configured URL
- Increase `DBAASP_TIMEOUT` if requests are timing out

### Data Validation Errors
- Verify CSV files have required columns (`Peptide ID` or `ID`, `N TERMINUS`)
- Check for proper UTF-8 encoding
- Ensure peptide IDs are numeric or in format `PREFIX_123`

### Enable Debug Logging
```bash
export LOG_LEVEL="DEBUG"
python main_api.py
```

## Example Workflow

```bash
# 1. Prepare your input file: data/input/peptides_C16.csv

# - data/output/lipophilicity_C16.csv
# 2. Run the analysis pipeline
python main_api.py

# 3. Check the generated output files
# - data/output/physchem_C16.csv
# - data/output/activity_C16.csv
# - data/output/activity_normalized_C16.csv
# - data/output/unified_results_C16.csv

# 4. Optionally run analysis scripts
python analysis/CorrelationMatrix.py
python analysis/HistSpecies.py
```

## Contributing

When modifying the pipeline:

1. Add proper logging statements using the configured logger
2. Raise appropriate custom exceptions from `src.core.exceptions`
3. Update configuration in `config.py` for any new parameters
4. Test with multiple N-terminal modifications

## License

[Add your license information here]

## Support

For issues or questions:
- Check the logs in `logs/pipeline.log`
- Enable DEBUG logging for detailed information
- Review the source code documentation

## Lipophilicity Calculation

The pipeline calculates partition and distribution coefficients for peptides:

- **logP (Partition Coefficient)**: Measures lipophilicity of the peptide based on its SMILES representation
- **logD (Distribution Coefficient)**: pH-dependent lipophilicity accounting for ionization state
- **SMILES**: Simplified Molecular Input Line Entry System representation of the peptide structure

For peptides with N-terminal modifications (C12, C16):
- The base peptide sequence is converted to SMILES using RDKit
- The N-terminal modification (fatty acid chain) is represented as a SMILES fragment
- Both are combined to generate the complete SMILES for the modified peptide
- logP and logD are calculated from the combined SMILES using molecular descriptor libraries

**N-Terminal Modifications Supported:**
- **C12**: Dodecanoyl (12-carbon fatty acid) - SMILES: `CCCCCCCCCCCC(=O)`
- **C16**: Hexadecanoyl (16-carbon fatty acid) - SMILES: `CCCCCCCCCCCCCCCC(=O)`

The output `lipophilicity_{NTERMINUS}.csv` contains:
- Peptide ID, sequence, N-terminus modification
- Generated SMILES for the modified peptide
- Calculated logP and logD values
