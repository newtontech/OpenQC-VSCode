"""Standard error types for the OpenQC Python bridge protocol."""

class BridgeError(Exception):
    """Base error for bridge operations."""
    pass

class MissingDependencyError(BridgeError):
    """A required Python dependency is not installed."""
    def __init__(self, package: str, install_hint: str = ""):
        self.package = package
        self.install_hint = install_hint
        super().__init__(f"Missing dependency: {package}. {install_hint}".strip())

class ParseError(BridgeError):
    """Failed to parse input file."""
    pass

class BridgeProtocolError(BridgeError):
    """Invalid bridge protocol message."""
    pass
