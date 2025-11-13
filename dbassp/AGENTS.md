# AGENTS.md - Developer Guidelines

## Build/Test/Lint Commands
- **Main execution**: `python main_api.py` (runs physchem, activity, normalize_activity pipeline)
- **Individual modules**: `python physchem.py`, `python activity.py`, `python normalize_activity.py`
- **No formal test framework configured** - verify by running scripts and checking output files
- **No linting configuration found** - follow PEP 8 standards

## Code Style Guidelines

### Imports
- Standard library first, then third-party (requests, pandas, matplotlib, seaborn)
- Use relative imports for local modules: `import physchem`, `import activity`

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
- Use pandas for CSV operations and data analysis
- Encode CSV files with "utf-8-sig" for proper BOM handling
- Use DictReader/DictWriter for structured CSV processing
- Default to empty strings for missing values rather than None