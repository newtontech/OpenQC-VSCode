#!/usr/bin/env python3
"""
OpenQC Format Converter - Python backend for quantum chemistry file format conversion.

This module uses dpdata library to convert between various quantum chemistry file formats.
Supported formats: VASP, Gaussian, ORCA, XYZ, PDB, CIF, and more.
"""

import sys
import json
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional

logger = logging.getLogger(__name__)


def convert_format(
    input_file: str,
    output_file: str,
    from_format: Optional[str] = None,
    to_format: Optional[str] = None,
    preserve_metadata: bool = True,
) -> Dict[str, Any]:
    """
    Convert quantum chemistry file between different formats.

    Args:
        input_file: Path to input file
        output_file: Path to output file
        from_format: Source format (auto-detected if None)
        to_format: Target format (inferred from output extension if None)
        preserve_metadata: Whether to preserve calculation metadata

    Returns:
        Dictionary with conversion result and metadata
    """
    try:
        import dpdata
    except ImportError:
        return {
            "success": False,
            "error": "dpdata library not found. Install with: pip install dpdata",
            "error_type": "ImportError",
        }

    input_path = Path(input_file)
    output_path = Path(output_file)

    if not input_path.exists():
        return {
            "success": False,
            "error": f"Input file not found: {input_file}",
            "error_type": "FileNotFoundError",
        }

    try:
        # Auto-detect input format from file extension if not specified
        if from_format is None:
            from_format = _detect_format(input_path)

        # Determine output format from extension if not specified
        if to_format is None:
            to_format = _detect_format(output_path)

        # Load system using dpdata
        system = _load_system(input_path, from_format)

        if system is None:
            return {
                "success": False,
                "error": f"Failed to load input file with format: {from_format}",
                "error_type": "LoadError",
            }

        # Extract metadata if requested
        metadata = {}
        if preserve_metadata:
            metadata = _extract_metadata(system, from_format)

        # Convert and save to output format
        _save_system(system, output_path, to_format)

        return {
            "success": True,
            "input_format": from_format,
            "output_format": to_format,
            "atoms_count": system.get_natoms(),
            "frames_count": len(system),
            "metadata": metadata,
            "output_file": str(output_path.absolute()),
        }

    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "error_type": type(e).__name__,
        }


def batch_convert(
    input_files: List[str],
    output_dir: str,
    to_format: str,
    from_format: Optional[str] = None,
    preserve_metadata: bool = True,
) -> Dict[str, Any]:
    """
    Batch convert multiple files to a target format.

    Args:
        input_files: List of input file paths
        output_dir: Directory for output files
        to_format: Target format
        from_format: Source format (auto-detected if None)
        preserve_metadata: Whether to preserve metadata

    Returns:
        Dictionary with batch conversion results
    """
    results = []
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for input_file in input_files:
        input_path = Path(input_file)
        output_file = output_path / f"{input_path.stem}.{_get_extension(to_format)}"

        result = convert_format(
            input_file=str(input_path),
            output_file=str(output_file),
            from_format=from_format,
            to_format=to_format,
            preserve_metadata=preserve_metadata,
        )

        results.append({
            "input_file": input_file,
            "output_file": str(output_file),
            "result": result,
        })

    successful = sum(1 for r in results if r["result"].get("success"))
    total = len(results)

    return {
        "success": successful == total,
        "total": total,
        "successful": successful,
        "failed": total - successful,
        "results": results,
    }


def supported_formats() -> Dict[str, List[str]]:
    """
    Get list of supported formats for conversion.

    Returns:
        Dictionary with input and output format lists
    """
    try:
        import dpdata
        # dpdata supports various formats
        return {
            "input": [
                "vasp", "poscar", "contcar",
                "gaussian", "gjf", "com",
                "orca", "inp",
                "xyz",
                "pdb",
                "cif",
                "cp2k",
                "qe", "quantum_espresso", "pw_in",
            ],
            "output": [
                "vasp", "poscar",
                "gaussian", "gjf", "com",
                "orca", "inp",
                "xyz",
                "pdb",
                "cif",
            ],
        }
    except ImportError:
        return {"input": [], "output": []}


def _detect_format(file_path: Path) -> str:
    """Auto-detect format from file extension and name."""
    filename = file_path.name.lower()
    ext = file_path.suffix.lower()

    # Special VASP files
    if filename == "poscar" or filename == "contcar":
        return "vasp/poscar"
    if filename == "incar":
        return "vasp/incar"
    if filename == "kpoints":
        return "vasp/kpoints"
    if filename == "potcar":
        return "vasp/potcar"

    # Extension-based detection
    format_map = {
        ".xyz": "xyz",
        ".pdb": "pdb",
        ".cif": "cif",
        ".gjf": "gaussian",
        ".com": "gaussian",
        ".inp": "orca",
    }

    return format_map.get(ext, "vasp/poscar")


def _get_extension(format_name: str) -> str:
    """Get file extension for a format."""
    ext_map = {
        "vasp": "poscar",
        "poscar": "poscar",
        "gaussian": "gjf",
        "gjf": "gjf",
        "com": "com",
        "orca": "inp",
        "inp": "inp",
        "xyz": "xyz",
        "pdb": "pdb",
        "cif": "cif",
    }
    return ext_map.get(format_name.lower(), "out")


def _load_system(file_path: Path, format_name: str):
    """Load system using dpdata."""
    import dpdata

    format_lower = format_name.lower()

    if "vasp" in format_lower or "poscar" in format_lower:
        return dpdata.System(str(file_path), fmt='vasp/poscar')
    elif "gaussian" in format_lower:
        return dpdata.System(str(file_path), fmt='gaussian')
    elif "orca" in format_lower:
        return dpdata.System(str(file_path), fmt='orca')
    elif "xyz" in format_lower:
        return dpdata.System(str(file_path), fmt='xyz')
    elif "pdb" in format_lower:
        return dpdata.System(str(file_path), fmt='pdb')
    elif "cif" in format_lower:
        return dpdata.System(str(file_path), fmt='cif')
    elif "cp2k" in format_lower:
        return dpdata.System(str(file_path), fmt='cp2k')
    else:
        # Try auto-detection
        return dpdata.System(str(file_path))


def _save_system(system, file_path: Path, format_name: str) -> None:
    """Save system using dpdata."""
    format_lower = format_name.lower()

    if "vasp" in format_lower or "poscar" in format_lower:
        system.to_vasp_poscar(str(file_path))
    elif "gaussian" in format_lower:
        system.to_gaussian_log(str(file_path))
    elif "xyz" in format_lower:
        system.to_xyz(str(file_path))
    elif "pdb" in format_lower:
        system.to_pdb(str(file_path))
    elif "cif" in format_lower:
        system.to_cif(str(file_path))
    else:
        # Default to POSCAR
        system.to_vasp_poscar(str(file_path))


def _extract_metadata(system, format_name: str) -> Dict[str, Any]:
    """Extract metadata from the system."""
    metadata = {
        "natoms": int(system.get_natoms()),
        "nbonds": int(getattr(system, "nbonds", 0)) if hasattr(system, "nbonds") else None,
    }

    # Add element information
    if hasattr(system, "data"):
        data = system.data
        if "atom_names" in data:
            metadata["elements"] = list(data["atom_names"])
        if "atom_numbs" in data:
            # Convert numpy int64 to regular int for JSON serialization
            metadata["atom_counts"] = [int(x) for x in data["atom_numbs"]]

    # Format-specific metadata
    if "vasp" in format_name.lower():
        metadata["original_format"] = "VASP"
    elif "gaussian" in format_name.lower():
        metadata["original_format"] = "Gaussian"
    elif "orca" in format_name.lower():
        metadata["original_format"] = "ORCA"

    return metadata


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Convert quantum chemistry file formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert VASP POSCAR to XYZ
  python format_converter.py POSCAR output.xyz

  # Convert Gaussian to ORCA
  python format_converter.py input.gjf output.inp --from gaussian --to orca

  # Batch convert multiple files
  python format_converter.py --batch file1.xyz file2.xyz --to vasp --output-dir ./converted

  # List supported formats
  python format_converter.py --list-formats
        """
    )

    parser.add_argument("input", nargs="?", help="Input file path")
    parser.add_argument("output", nargs="?", help="Output file path")
    parser.add_argument("--from", dest="from_format", help="Input format")
    parser.add_argument("--to", dest="to_format", help="Output format")
    parser.add_argument("--batch", nargs="+", help="Batch convert multiple files")
    parser.add_argument("--output-dir", help="Output directory for batch conversion")
    parser.add_argument("--no-metadata", action="store_true", help="Don't preserve metadata")
    parser.add_argument("--list-formats", action="store_true", help="List supported formats")
    parser.add_argument("--json", action="store_true", help="Output result as JSON")

    args = parser.parse_args()

    if args.list_formats:
        formats = supported_formats()
        print("Supported input formats:", ", ".join(formats["input"]))
        print("Supported output formats:", ", ".join(formats["output"]))
        return 0

    if args.batch:
        if not args.to_format:
            logger.error("--to format required for batch conversion")
            print("Error: --to format required for batch conversion", file=sys.stderr)
            return 1

        result = batch_convert(
            input_files=args.batch,
            output_dir=args.output_dir or "./converted",
            to_format=args.to_format,
            from_format=args.from_format,
            preserve_metadata=not args.no_metadata,
        )
    else:
        if not args.input or not args.output:
            parser.print_help()
            return 1

        result = convert_format(
            input_file=args.input,
            output_file=args.output,
            from_format=args.from_format,
            to_format=args.to_format,
            preserve_metadata=not args.no_metadata,
        )

    if args.json or args.batch:
        print(json.dumps(result, indent=2))
    else:
        if result.get("success"):
            logger.info(f"Successfully converted to {result['output_format']} format")
            print(f"Successfully converted to {result['output_format']} format")
            print(f"Output: {result.get('output_file')}")
            print(f"Atoms: {result.get('atoms_count')}")
        else:
            logger.error(f"Conversion failed: {result.get('error')}")
            print(f"Error: {result.get('error')}", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
