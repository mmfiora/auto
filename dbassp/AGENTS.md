# AGENTS.md - Developer Guidelines

## Project Overview
DBAASP peptide data pipeline: fetches raw physicochemical and activity data from the
DBAASP API, normalizes concentrations, and writes output CSVs.

**Output files live in Google Drive** via a directory junction:
- `data/output/` → `G:\My Drive\Doctorado\Automatizacion\auto_sqlite\data\`
- `presentations/` → `G:\My Drive\Doctorado\Automatizacion\auto_sqlite\presentations\`

Run `python scripts/setup_symlinks.py` (as Admin) once to create these junctions.

## Build/Test/Lint Commands
- **Full pipeline**: `python main_api.py` (interactive) or `python main_api.py --nterminus C16`
- **All N-termini**: `python scripts/collection/fetch_all_to_csv.py`
- **Database build**: `python scripts/database/create_sqlite_db.py`
- **No formal test framework** — verify by running scripts and checking output CSVs
- **No linting configuration** — follow PEP 8 standards

## Pipeline Steps (in order)
1. `physchem.run()` — fetch raw physicochemical properties from DBAASP API
2. `activity.run()` — fetch raw activity data + resolve full citations (`reference_citation` column)
3. `normalize_activity.run()` — parse concentrations → µg/mL, classify activity, split species/strain
4. `intrinsic_properties.run()` — combine physchem + derived scalars (long_tail, MW, total_charge) + min-list
5. `unified_results.run()` — join normalized activity with physchem properties
6. `activity_summary.run()` — summary statistics per peptide
7. `generate_peptide_list.run()` — generate peptide list outputs

> **Note**: logP / logD are NOT calculated here. Use a separate project for lipophilicity.

## Reference Citation Column
`activity.py` adds a `reference_citation` column to the activity CSV. Each activity row in the
DBAASP API has a `reference` field (a 1-based index string, e.g. `"1"`) that points into the
peptide's top-level `articles` list. The citation is formatted as:
`"Authors | Title | Journal Year;Vol:Pages | PMID:xxx"`

## Code Style Guidelines

### Imports
- Standard library first, then third-party (requests, pandas, matplotlib, seaborn)
- Use package imports for local modules: `from src.core import common`

### Formatting & Structure
- Use double quotes for strings consistently
- 4-space indentation (no tabs)
- Keep functions focused and well-documented with docstrings
- Use type hints where helpful: `def fetch(pid: int):`, `def calc_mw(seq: str, nterm: str | None, cterm: str | None) -> float:`

### Naming Conventions
- **Files**: snake_case (e.g., `main_api.py`, `normalize_activity.py`)
- **Functions**: snake_case (e.g., `load_ids()`, `calc_mw()`)
- **Variables**: snake_case (e.g., `api_url`, `corr_matrix`)
- **Constants**: UPPER_CASE (e.g., `API_URL`, `HEADERS`, `AA_MASS`)

### Error Handling
- Use try/except blocks for API calls and file operations
- Print informative error messages: `print(f" fail ({e})")`
- Use `resp.raise_for_status()` for HTTP requests
- Gracefully handle missing data with empty string defaults

### Data Processing
- Encode CSV files with `"utf-8-sig"` for proper BOM handling
- Use `DictReader`/`DictWriter` for structured CSV processing
- Default to empty strings for missing values rather than None