"""
Custom exceptions for the Yleaf package.
"""

class YleafError(Exception):
    """Base exception for all Yleaf errors."""
    pass

class ConfigurationError(YleafError):
    """Raised when there is an issue with the configuration."""
    pass

class ExternalCommandError(YleafError):
    """Raised when an external command (subprocess) fails."""
    def __init__(self, command: str, return_code: int, stderr: str):
        self.command = command
        self.return_code = return_code
        self.stderr = stderr
        super().__init__(f"Command '{command}' failed with return code {return_code}. Error: {stderr}")

class InputError(YleafError):
    """Raised when input data is invalid."""
    pass
