"""
ASE-based format converter.
"""

from pathlib import Path
from typing import Optional, Union, List


class ASEConverter:
    """
    Format converter using ASE (Atomic Simulation Environment).
    
    Supports all ASE formats including:
    - XYZ
    - CIF
    - PDB
    - VASP
    - Quantum ESPRESSO
    - And many more
    """
    
    # Format mapping
    FORMAT_MAP = {
        '.xyz': 'xyz',
        '.cif': 'cif',
        '.pdb': 'pdb',
        '.vasp': 'vasp',
        'poscar': 'vasp',
        'contcar': 'vasp',
        '.pw': 'espresso',
        '.scf': 'espresso',
        '.mol': 'mol',
        '.mol2': 'mol2',
        '.json': 'json',
    }
    
    def __init__(self):
        self._ase_io = None
        self._check_ase()
    
    def _check_ase(self) -> None:
        """Check if ASE is available."""
        try:
            from ase import io
            self._ase_io = io
        except ImportError:
            raise ImportError(
                "ASE is not installed. Install with:\n"
                "  pip install ase"
            )
    
    def convert(
        self,
        input_file: Union[str, Path],
        output_file: Union[str, Path],
        output_format: Optional[str] = None,
        **kwargs
    ) -> Path:
        """
        Convert input file to output format.
        
        Args:
            input_file: Path to input file
            output_file: Path to output file
            output_format: Target format (auto-detected if not provided)
            **kwargs: Additional arguments passed to ASE
            
        Returns:
            Path to the converted file
        """
        input_path = Path(input_file)
        output_path = Path(output_file)
        
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")
        
        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Detect input format
        input_format = self._detect_format(input_path)
        
        # Detect output format
        if output_format is None:
            output_format = self._detect_format(output_path, for_output=True)
        
        # Read input
        atoms = self._ase_io.read(str(input_path), format=input_format, **kwargs)
        
        # Write output
        self._ase_io.write(str(output_path), atoms, format=output_format)
        
        return output_path
    
    def _detect_format(
        self, 
        filepath: Path, 
        for_output: bool = False
    ) -> str:
        """Detect file format from extension/name."""
        suffix = filepath.suffix.lower()
        name = filepath.name.lower()
        
        # Check by name first
        if name in self.FORMAT_MAP:
            return self.FORMAT_MAP[name]
        
        # Check by extension
        if suffix in self.FORMAT_MAP:
            return self.FORMAT_MAP[suffix]
        
        # Special cases
        if name in ['poscar', 'contcar']:
            return 'vasp'
        
        # Default to extension name
        return suffix.lstrip('.')
    
    def get_supported_formats(self) -> List[str]:
        """Get list of supported formats."""
        return list(set(self.FORMAT_MAP.values()))
    
    def read_structure(self, filepath: Union[str, Path]):
        """Read structure from file and return ASE Atoms object."""
        from ase import io
        return io.read(str(filepath))
    
    def visualize(self, filepath: Union[str, Path], viewer: str = 'ase'):
        """Visualize structure using ASE viewer."""
        atoms = self.read_structure(filepath)
        atoms.edit()


# Singleton instance
ase_converter = ASEConverter()
