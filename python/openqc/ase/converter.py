"""
ASE Converter - Format Conversion Core

Provides conversion between quantum chemistry file formats and ASE Atoms objects.
Supports VASP, CP2K, QE, Gaussian, ORCA, NWChem, GAMESS, and LAMMPS formats.
"""

from typing import Dict, Any, Optional, List, Union
from pathlib import Path
import json
import sys
from ase import Atoms
from ase.io import read, write
from ase.io.vasp import read_vasp, write_vasp
from ase.io.cp2k import read_cp2k_input
import numpy as np


class ConversionResult:
    """Result of a conversion operation."""
    
    def __init__(
        self,
        success: bool,
        data: Optional[Union[Atoms, str, Dict]] = None,
        error: Optional[str] = None,
        warnings: List[str] = None,
        metadata: Dict[str, Any] = None
    ):
        self.success = success
        self.data = data
        self.error = error
        self.warnings = warnings or []
        self.metadata = metadata or {}
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary for JSON serialization."""
        result = {
            'success': self.success,
            'error': self.error,
            'warnings': self.warnings,
            'metadata': self.metadata
        }
        
        if isinstance(self.data, Atoms):
            result['data_type'] = 'atoms'
            result['atoms'] = self._atoms_to_dict(self.data)
        elif isinstance(self.data, str):
            result['data_type'] = 'text'
            result['content'] = self.data
        elif isinstance(self.data, dict):
            result['data_type'] = 'dict'
            result['data'] = self.data
        
        return result
    
    def _atoms_to_dict(self, atoms: Atoms) -> Dict[str, Any]:
        """Convert ASE Atoms to dictionary."""
        data = {
            'chemical_symbols': atoms.get_chemical_symbols(),
            'positions': atoms.get_positions().tolist(),
            'cell': atoms.get_cell().tolist() if atoms.pbc.any() else None,
            'pbc': atoms.pbc.tolist(),
            'info': dict(atoms.info) if atoms.info else {}
        }
        
        if atoms.has('masses'):
            data['masses'] = atoms.get_masses().tolist()
        
        if atoms.constraints:
            data['constraints'] = str(atoms.constraints)
        
        return data


class ASEConverter:
    """
    ASE-based format converter for quantum chemistry files.
    
    Supports bidirectional conversion between:
    - VASP (POSCAR/CONTCAR)
    - CP2K (.inp)
    - Quantum ESPRESSO (.in, .pw.in)
    - Gaussian (.com, .gjf)
    - ORCA (.inp)
    - NWChem (.nw, .nwinp)
    - GAMESS (.inp)
    - LAMMPS (.lmp, .data)
    - XYZ (.xyz)
    - PDB (.pdb)
    - CIF (.cif)
    """
    
    def __init__(self):
        self.supported_formats = {
            'vasp': {'extensions': ['POSCAR', 'CONTCAR'], 'ase_format': 'vasp'},
            'cp2k': {'extensions': ['.inp'], 'ase_format': 'cp2k'},
            'qe': {'extensions': ['.in', '.pw.in'], 'ase_format': 'espresso'},
            'gaussian': {'extensions': ['.com', '.gjf'], 'ase_format': 'gaussian'},
            'orca': {'extensions': ['.inp'], 'ase_format': 'orca'},
            'nwchem': {'extensions': ['.nw', '.nwinp'], 'ase_format': 'nwchem'},
            'gamess': {'extensions': ['.inp'], 'ase_format': 'gamess'},
            'lammps': {'extensions': ['.lmp', '.data'], 'ase_format': 'lammps-data'},
            'xyz': {'extensions': ['.xyz'], 'ase_format': 'xyz'},
            'pdb': {'extensions': ['.pdb'], 'ase_format': 'pdb'},
            'cif': {'extensions': ['.cif'], 'ase_format': 'cif'}
        }
    
    def read_to_atoms(
        self,
        filepath: str,
        format_hint: Optional[str] = None,
        index: int = -1
    ) -> ConversionResult:
        """
        Read file to ASE Atoms object.
        
        Parameters
        ----------
        filepath : str
            Path to input file
        format_hint : Optional[str]
            Format hint (vasp, cp2k, qe, gaussian, etc.)
        index : int
            Index of structure to read (-1 for last, ':' for all)
        
        Returns
        -------
        ConversionResult
            Result with Atoms object or error
        """
        try:
            path = Path(filepath)
            if not path.exists():
                return ConversionResult(
                    success=False,
                    error=f"File not found: {filepath}"
                )
            
            # Determine format
            ase_format = self._detect_format(filepath, format_hint)
            if not ase_format:
                return ConversionResult(
                    success=False,
                    error=f"Could not determine format for: {filepath}"
                )
            
            # Read file
            atoms = read(filepath, format=ase_format, index=index)
            
            if isinstance(atoms, list):
                # Multiple structures - return the last one or specified index
                if len(atoms) == 0:
                    return ConversionResult(
                        success=False,
                        error="No structures found in file"
                    )
                atoms = atoms[-1] if index == -1 else atoms[index]
            
            return ConversionResult(
                success=True,
                data=atoms,
                metadata={
                    'source_file': filepath,
                    'source_format': ase_format,
                    'natoms': len(atoms),
                    'formula': atoms.get_chemical_formula()
                }
            )
            
        except Exception as e:
            return ConversionResult(
                success=False,
                error=f"Failed to read file: {str(e)}"
            )
    
    def write_from_atoms(
        self,
        atoms: Atoms,
        output_path: str,
        format_hint: Optional[str] = None,
        **kwargs
    ) -> ConversionResult:
        """
        Write ASE Atoms object to file.
        
        Parameters
        ----------
        atoms : Atoms
            ASE Atoms object to write
        output_path : str
            Output file path
        format_hint : Optional[str]
            Format hint (vasp, cp2k, qe, gaussian, etc.)
        **kwargs
            Additional format-specific parameters
        
        Returns
        -------
        ConversionResult
            Result with output path or error
        """
        try:
            if not isinstance(atoms, Atoms):
                return ConversionResult(
                    success=False,
                    error="Input is not an ASE Atoms object"
                )
            
            # Determine format
            ase_format = self._detect_format(output_path, format_hint)
            if not ase_format:
                return ConversionResult(
                    success=False,
                    error=f"Could not determine format for: {output_path}"
                )
            
            # Write file
            write(output_path, atoms, format=ase_format, **kwargs)
            
            return ConversionResult(
                success=True,
                data=output_path,
                metadata={
                    'output_file': output_path,
                    'output_format': ase_format,
                    'natoms': len(atoms)
                }
            )
            
        except Exception as e:
            return ConversionResult(
                success=False,
                error=f"Failed to write file: {str(e)}"
            )
    
    def convert_format(
        self,
        input_path: str,
        output_path: str,
        input_format: Optional[str] = None,
        output_format: Optional[str] = None,
        **kwargs
    ) -> ConversionResult:
        """
        Convert between two formats.
        
        Parameters
        ----------
        input_path : str
            Input file path
        output_path : str
            Output file path
        input_format : Optional[str]
            Input format hint
        output_format : Optional[str]
            Output format hint
        **kwargs
            Additional format-specific parameters
        
        Returns
        -------
        ConversionResult
            Result with conversion metadata
        """
        # Read input
        read_result = self.read_to_atoms(input_path, format_hint=input_format)
        if not read_result.success:
            return read_result
        
        atoms = read_result.data
        
        # Write output
        write_result = self.write_from_atoms(
            atoms,
            output_path,
            format_hint=output_format,
            **kwargs
        )
        
        if not write_result.success:
            return write_result
        
        return ConversionResult(
            success=True,
            data=output_path,
            metadata={
                **read_result.metadata,
                **write_result.metadata,
                'conversion': f"{read_result.metadata['source_format']} -> {write_result.metadata['output_format']}"
            },
            warnings=read_result.warnings + write_result.warnings
        )
    
    def _detect_format(
        self,
        filepath: str,
        format_hint: Optional[str] = None
    ) -> Optional[str]:
        """
        Detect file format from path and hint.
        
        Parameters
        ----------
        filepath : str
            File path
        format_hint : Optional[str]
            User-provided format hint
        
        Returns
        -------
        Optional[str]
            ASE format identifier or None
        """
        # Use hint if provided
        if format_hint:
            hint_lower = format_hint.lower()
            if hint_lower in self.supported_formats:
                return self.supported_formats[hint_lower]['ase_format']
        
        # Detect from filename
        path = Path(filepath)
        filename = path.name.upper()
        extension = path.suffix.lower()
        
        # Special cases for VASP
        if filename in ['POSCAR', 'CONTCAR']:
            return 'vasp'
        
        # Check by extension
        for fmt_name, fmt_info in self.supported_formats.items():
            if extension in fmt_info['extensions']:
                return fmt_info['ase_format']
            if filename in fmt_info['extensions']:
                return fmt_info['ase_format']
        
        # Try common formats
        common_extensions = {
            '.xyz': 'xyz',
            '.pdb': 'pdb',
            '.cif': 'cif',
            '.com': 'gaussian',
            '.gjf': 'gaussian',
            '.inp': None,  # Ambiguous - could be ORCA, CP2K, etc.
        }
        
        return common_extensions.get(extension)


def main():
    """CLI interface for ASE converter."""
    if len(sys.argv) < 2:
        print("Usage: python converter.py <command> [args]")
        print("Commands:")
        print("  read <filepath> [format]")
        print("  write <atoms_json> <output_path> [format]")
        print("  convert <input_path> <output_path> [input_format] [output_format]")
        sys.exit(1)
    
    command = sys.argv[1]
    converter = ASEConverter()
    
    if command == 'read':
        if len(sys.argv) < 3:
            print("Error: filepath required")
            sys.exit(1)
        
        filepath = sys.argv[2]
        format_hint = sys.argv[3] if len(sys.argv) > 3 else None
        
        result = converter.read_to_atoms(filepath, format_hint)
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == 'convert':
        if len(sys.argv) < 4:
            print("Error: input_path and output_path required")
            sys.exit(1)
        
        input_path = sys.argv[2]
        output_path = sys.argv[3]
        input_format = sys.argv[4] if len(sys.argv) > 4 else None
        output_format = sys.argv[5] if len(sys.argv) > 5 else None
        
        result = converter.convert_format(
            input_path,
            output_path,
            input_format,
            output_format
        )
        print(json.dumps(result.to_dict(), indent=2))
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)


if __name__ == '__main__':
    main()
