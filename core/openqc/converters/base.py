"""
Base converter functionality.
"""

from pathlib import Path
from typing import Optional, Union
from openqc.parsers.base import parse_file, ParsedInput


def convert_format(
    input_file: Union[str, Path],
    output_file: Union[str, Path],
    output_format: Optional[str] = None,
) -> Path:
    """
    Convert between quantum chemistry file formats.
    
    Args:
        input_file: Path to input file
        output_file: Path to output file
        output_format: Target format (auto-detected from output_file if not provided)
        
    Returns:
        Path to the converted file
    """
    input_path = Path(input_file)
    output_path = Path(output_file)
    
    # Auto-detect output format from extension
    if output_format is None:
        output_format = output_path.suffix.lower().lstrip('.')
    
    # Try dpdata first (most comprehensive)
    try:
        from openqc.converters.dpdata_converter import dpdata_converter
        return dpdata_converter.convert(input_path, output_path, output_format)
    except ImportError:
        pass
    
    # Fall back to ASE
    try:
        from openqc.converters.ase_converter import ase_converter
        return ase_converter.convert(input_path, output_path, output_format)
    except ImportError:
        pass
    
    raise RuntimeError(
        "No converter available. Please install dpdata or ase:\n"
        "  pip install dpdata\n"
        "  pip install ase"
    )
