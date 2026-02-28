"""
Base parser class for quantum chemistry input files.
"""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
from pydantic import BaseModel


class Atom(BaseModel):
    """Represents an atom in the structure."""
    element: str
    x: float
    y: float
    z: float
    index: Optional[int] = None


class Structure(BaseModel):
    """Represents a molecular or crystal structure."""
    atoms: List[Atom]
    cell: Optional[List[List[float]]] = None
    title: Optional[str] = None
    charge: int = 0
    multiplicity: int = 1
    
    @property
    def num_atoms(self) -> int:
        return len(self.atoms)
    
    @property
    def elements(self) -> List[str]:
        return [atom.element for atom in self.atoms]
    
    def get_positions(self) -> List[List[float]]:
        return [[atom.x, atom.y, atom.z] for atom in self.atoms]


class CalculationParams(BaseModel):
    """Calculation parameters."""
    method: Optional[str] = None
    basis: Optional[str] = None
    job_type: Optional[str] = None
    functional: Optional[str] = None
    additional_keywords: Dict[str, Any] = {}


class ParsedInput(BaseModel):
    """Result of parsing a quantum chemistry input file."""
    structure: Optional[Structure] = None
    parameters: Optional[CalculationParams] = None
    raw_content: str = ""
    file_type: str = ""
    errors: List[str] = []
    warnings: List[str] = []


class BaseParser(ABC):
    """Abstract base class for all parsers."""
    
    def __init__(self):
        self.errors: List[str] = []
        self.warnings: List[str] = []
    
    @abstractmethod
    def parse(self, content: str) -> ParsedInput:
        """Parse the content and return a ParsedInput object."""
        pass
    
    def parse_file(self, filepath: Union[str, Path]) -> ParsedInput:
        """Parse a file and return a ParsedInput object."""
        filepath = Path(filepath)
        
        if not filepath.exists():
            raise FileNotFoundError(f"File not found: {filepath}")
        
        content = filepath.read_text()
        result = self.parse(content)
        result.raw_content = content
        
        return result
    
    def add_error(self, message: str) -> None:
        """Add an error message."""
        self.errors.append(message)
    
    def add_warning(self, message: str) -> None:
        """Add a warning message."""
        self.warnings.append(message)
    
    def reset(self) -> None:
        """Reset error and warning lists."""
        self.errors = []
        self.warnings = []


def parse_file(filepath: Union[str, Path]) -> ParsedInput:
    """
    Auto-detect file type and parse.
    
    Args:
        filepath: Path to the input file
        
    Returns:
        ParsedInput object
    """
    filepath = Path(filepath)
    suffix = filepath.suffix.lower()
    name = filepath.name.lower()
    
    # Determine parser based on file extension/name
    if suffix in ['.gjf', '.com', '.gau']:
        from openqc.parsers.gaussian import GaussianParser
        parser = GaussianParser()
    elif suffix == '.vasp' or name in ['poscar', 'contcar']:
        from openqc.parsers.vasp import VASPParser
        parser = VASPParser()
    elif suffix in ['.pw', '.scf', '.nscf', '.relax', '.md']:
        from openqc.parsers.qe import QEParser
        parser = QEParser()
    elif suffix == '.inp':
        # Could be ORCA or CP2K, try ORCA first
        from openqc.parsers.orca import ORCAParser
        parser = ORCAParser()
    elif suffix == '.xyz':
        from openqc.parsers.xyz import XYZParser
        parser = XYZParser()
    else:
        raise ValueError(f"Unknown file type: {suffix}")
    
    return parser.parse_file(filepath)
