"""
dpdata-based format converter.
"""

from pathlib import Path
from typing import Optional, Union, List


class DpdataConverter:
    """
    Format converter using dpdata library.
    
    Supports 50+ formats including:
    - Gaussian (.gjf, .log)
    - VASP (POSCAR, CONTCAR)
    - Quantum ESPRESSO
    - ORCA
    - GROMACS
    - LAMMPS
    - CP2K
    - And many more
    """
    
    # Format mapping
    FORMAT_MAP = {
        '.gjf': 'gaussian/gjf',
        '.com': 'gaussian/gjf',
        '.log': 'gaussian/log',
        '.vasp': 'vasp/poscar',
        'poscar': 'vasp/poscar',
        'contcar': 'vasp/poscar',
        '.pw': 'pwscf',
        '.scf': 'pwscf',
        '.inp': 'orca',
        '.xyz': 'xyz',
        '.cif': 'cif',
        '.pdb': 'pdb',
        '.mol2': 'mol2',
        '.gro': 'gromacs/gro',
        '.lammpstrj': 'lammps/lmp',
        '.lmp': 'lammps/lmp',
    }
    
    def __init__(self):
        self._dpdata = None
        self._check_dpdata()
    
    def _check_dpdata(self) -> None:
        """Check if dpdata is available."""
        try:
            import dpdata
            self._dpdata = dpdata
        except ImportError:
            raise ImportError(
                "dpdata is not installed. Install with:\n"
                "  pip install dpdata"
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
            **kwargs: Additional arguments passed to dpdata
            
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
        data = self._dpdata.System(str(input_path), fmt=input_format, **kwargs)
        
        # Write output
        data.to(output_format, str(output_path))
        
        return output_path
    
    def _detect_format(
        self, 
        filepath: Path, 
        for_output: bool = False
    ) -> str:
        """Detect file format from extension/name."""
        suffix = filepath.suffix.lower()
        name = filepath.name.lower()
        
        # Check by name first (for POSCAR, CONTCAR, etc.)
        if name in self.FORMAT_MAP:
            return self.FORMAT_MAP[name]
        
        # Check by extension
        if suffix in self.FORMAT_MAP:
            return self.FORMAT_MAP[suffix]
        
        # Special cases
        if name in ['poscar', 'contcar']:
            return 'vasp/poscar'
        if name == 'incar':
            return 'vasp/incar'
        
        # Default to extension name
        return suffix.lstrip('.')
    
    def get_supported_formats(self) -> List[str]:
        """Get list of supported formats."""
        return list(set(self.FORMAT_MAP.values()))
    
    def batch_convert(
        self,
        input_files: List[Union[str, Path]],
        output_dir: Union[str, Path],
        output_format: str,
    ) -> List[Path]:
        """
        Convert multiple files to a single format.
        
        Args:
            input_files: List of input file paths
            output_dir: Output directory
            output_format: Target format
            
        Returns:
            List of converted file paths
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        converted = []
        for input_file in input_files:
            input_path = Path(input_file)
            output_file = output_path / f"{input_path.stem}.{output_format}"
            converted.append(self.convert(input_path, output_file, output_format))
        
        return converted


# Singleton instance
dpdata_converter = DpdataConverter()
