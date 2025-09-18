# exceptions.py
# Custom exceptions for DBAASP pipeline

class DBAASSPError(Exception):
    """Base exception for DBAASP pipeline errors."""
    pass

class APIError(DBAASSPError):
    """Exception raised for API-related errors."""
    def __init__(self, message: str, peptide_id: int | None = None, status_code: int | None = None):
        self.peptide_id = peptide_id
        self.status_code = status_code
        super().__init__(message)

class FileProcessingError(DBAASSPError):
    """Exception raised for file processing errors."""
    def __init__(self, message: str, filename: str | None = None):
        self.filename = filename
        super().__init__(message)

class DataValidationError(DBAASSPError):
    """Exception raised for data validation errors."""
    def __init__(self, message: str, field: str | None = None, value: str | None = None):
        self.field = field
        self.value = value
        super().__init__(message)

class ConfigurationError(DBAASSPError):
    """Exception raised for configuration errors."""
    pass