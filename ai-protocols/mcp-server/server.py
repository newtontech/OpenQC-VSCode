#!/usr/bin/env python3
"""
OpenQC MCP Server

A Model Context Protocol (MCP) server for Claude Code integration.
Provides tools for parsing, converting, visualizing, and modifying
quantum chemistry input files.
"""

import asyncio
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

# MCP SDK imports
try:
    from mcp.server import Server
    from mcp.server.stdio import stdio_server
    from mcp.types import Tool, TextContent
except ImportError:
    print("MCP SDK not found. Install with: pip install mcp", file=sys.stderr)
    sys.exit(1)

# Create server instance
server = Server("openqc-mcp-server")


@server.list_tools()
async def list_tools() -> List[Tool]:
    """List available MCP tools."""
    return [
        Tool(
            name="parse_qc_file",
            description="Parse a quantum chemistry input file and extract structure and parameters",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the QC input file"
                    }
                },
                "required": ["file_path"]
            }
        ),
        Tool(
            name="convert_format",
            description="Convert between quantum chemistry file formats",
            inputSchema={
                "type": "object",
                "properties": {
                    "input_file": {
                        "type": "string",
                        "description": "Path to input file"
                    },
                    "output_file": {
                        "type": "string",
                        "description": "Path to output file"
                    },
                    "output_format": {
                        "type": "string",
                        "description": "Target format (e.g., 'vasp', 'xyz', 'cif')"
                    }
                },
                "required": ["input_file", "output_file"]
            }
        ),
        Tool(
            name="visualize_structure",
            description="Generate 3D visualization of molecular structure",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the structure file"
                    },
                    "style": {
                        "type": "string",
                        "enum": ["stick", "sphere", "cartoon", "line"],
                        "default": "stick",
                        "description": "Visualization style"
                    },
                    "output_format": {
                        "type": "string",
                        "enum": ["html", "json"],
                        "default": "html",
                        "description": "Output format for visualization"
                    }
                },
                "required": ["file_path"]
            }
        ),
        Tool(
            name="modify_parameters",
            description="Modify calculation parameters in a QC input file",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the input file"
                    },
                    "modifications": {
                        "type": "object",
                        "description": "Parameters to modify (e.g., {'functional': 'B3LYP', 'basis': 'def2-TZVP'})"
                    },
                    "output_file": {
                        "type": "string",
                        "description": "Path to save modified file (optional)"
                    }
                },
                "required": ["file_path", "modifications"]
            }
        ),
        Tool(
            name="validate_input",
            description="Validate a quantum chemistry input file for errors",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the input file"
                    }
                },
                "required": ["file_path"]
            }
        ),
        Tool(
            name="get_molecule_info",
            description="Get basic information about a molecule (formula, mass, etc.)",
            inputSchema={
                "type": "object",
                "properties": {
                    "file_path": {
                        "type": "string",
                        "description": "Path to the structure file"
                    }
                },
                "required": ["file_path"]
            }
        ),
        Tool(
            name="generate_input_template",
            description="Generate a template input file for a specific QC software",
            inputSchema={
                "type": "object",
                "properties": {
                    "software": {
                        "type": "string",
                        "enum": ["gaussian", "vasp", "orca", "qe"],
                        "description": "Target QC software"
                    },
                    "job_type": {
                        "type": "string",
                        "enum": ["sp", "opt", "freq", "md"],
                        "default": "sp",
                        "description": "Type of calculation"
                    },
                    "functional": {
                        "type": "string",
                        "default": "B3LYP",
                        "description": "DFT functional"
                    },
                    "basis": {
                        "type": "string",
                        "default": "def2-SVP",
                        "description": "Basis set"
                    }
                },
                "required": ["software"]
            }
        ),
    ]


@server.call_tool()
async def call_tool(name: str, arguments: Dict[str, Any]) -> List[TextContent]:
    """Handle tool calls."""
    
    try:
        if name == "parse_qc_file":
            return await handle_parse_qc_file(arguments)
        elif name == "convert_format":
            return await handle_convert_format(arguments)
        elif name == "visualize_structure":
            return await handle_visualize_structure(arguments)
        elif name == "modify_parameters":
            return await handle_modify_parameters(arguments)
        elif name == "validate_input":
            return await handle_validate_input(arguments)
        elif name == "get_molecule_info":
            return await handle_get_molecule_info(arguments)
        elif name == "generate_input_template":
            return await handle_generate_input_template(arguments)
        else:
            return [TextContent(type="text", text=f"Unknown tool: {name}")]
    except Exception as e:
        return [TextContent(type="text", text=f"Error: {str(e)}")]


async def handle_parse_qc_file(args: Dict[str, Any]) -> List[TextContent]:
    """Parse a QC input file."""
    file_path = Path(args["file_path"])
    
    if not file_path.exists():
        return [TextContent(type="text", text=f"File not found: {file_path}")]
    
    # Import and use the parser
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / "core"))
    from openqc.parsers import parse_file
    
    result = parse_file(file_path)
    
    output = {
        "file_type": result.file_type,
        "errors": result.errors,
        "warnings": result.warnings,
    }
    
    if result.structure:
        output["structure"] = {
            "num_atoms": result.structure.num_atoms,
            "elements": result.structure.elements,
            "charge": result.structure.charge,
            "multiplicity": result.structure.multiplicity,
            "title": result.structure.title,
        }
    
    if result.parameters:
        output["parameters"] = {
            "method": result.parameters.method,
            "basis": result.parameters.basis,
            "functional": result.parameters.functional,
            "job_type": result.parameters.job_type,
        }
    
    return [TextContent(type="text", text=json.dumps(output, indent=2))]


async def handle_convert_format(args: Dict[str, Any]) -> List[TextContent]:
    """Convert between formats."""
    input_file = Path(args["input_file"])
    output_file = Path(args["output_file"])
    output_format = args.get("output_format")
    
    if not input_file.exists():
        return [TextContent(type="text", text=f"Input file not found: {input_file}")]
    
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / "core"))
    from openqc.converters import convert_format
    
    result_path = convert_format(input_file, output_file, output_format)
    
    return [TextContent(
        type="text", 
        text=f"Successfully converted to: {result_path}"
    )]


async def handle_visualize_structure(args: Dict[str, Any]) -> List[TextContent]:
    """Generate visualization."""
    file_path = Path(args["file_path"])
    style = args.get("style", "stick")
    output_format = args.get("output_format", "html")
    
    if not file_path.exists():
        return [TextContent(type="text", text=f"File not found: {file_path}")]
    
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / "core"))
    from openqc.parsers import parse_file
    from openqc.visualizers import StructureVisualizer, export_visualization
    
    result = parse_file(file_path)
    
    if not result.structure:
        return [TextContent(type="text", text="No structure found in file")]
    
    viz = StructureVisualizer('3dmol')
    
    if output_format == 'json':
        data = viz.to_json(result.structure)
        return [TextContent(type="text", text=json.dumps(data, indent=2))]
    else:
        html = viz.visualize(result.structure, style=style)
        return [TextContent(type="text", text=html)]


async def handle_modify_parameters(args: Dict[str, Any]) -> List[TextContent]:
    """Modify input file parameters."""
    file_path = Path(args["file_path"])
    modifications = args["modifications"]
    output_file = args.get("output_file")
    
    if not file_path.exists():
        return [TextContent(type="text", text=f"File not found: {file_path}")]
    
    # Read file
    content = file_path.read_text()
    
    # Apply modifications (simplified implementation)
    # A full implementation would parse and regenerate the file
    modified = content
    
    for key, value in modifications.items():
        # Try to replace common patterns
        import re
        
        # Functional replacement
        if key.lower() == 'functional':
            for func in ['B3LYP', 'PBE0', 'M06-2X', 'wB97X-D', 'CAM-B3LYP']:
                modified = re.sub(
                    rf'\b{func}\b', 
                    str(value), 
                    modified, 
                    flags=re.IGNORECASE
                )
        
        # Basis set replacement
        if key.lower() == 'basis':
            for basis in ['6-31G(d)', 'def2-SVP', 'def2-TZVP', 'cc-pVDZ']:
                modified = re.sub(
                    rf'\b{re.escape(basis)}\b',
                    str(value),
                    modified,
                    flags=re.IGNORECASE
                )
    
    # Save
    if output_file:
        Path(output_file).write_text(modified)
        return [TextContent(type="text", text=f"Modified file saved to: {output_file}")]
    else:
        return [TextContent(type="text", text=modified)]


async def handle_validate_input(args: Dict[str, Any]) -> List[TextContent]:
    """Validate input file."""
    file_path = Path(args["file_path"])
    
    if not file_path.exists():
        return [TextContent(type="text", text=f"File not found: {file_path}")]
    
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / "core"))
    from openqc.parsers import parse_file
    
    result = parse_file(file_path)
    
    validation_result = {
        "valid": len(result.errors) == 0,
        "errors": result.errors,
        "warnings": result.warnings,
        "file_type": result.file_type,
    }
    
    return [TextContent(type="text", text=json.dumps(validation_result, indent=2))]


async def handle_get_molecule_info(args: Dict[str, Any]) -> List[TextContent]:
    """Get molecule information."""
    file_path = Path(args["file_path"])
    
    if not file_path.exists():
        return [TextContent(type="text", text=f"File not found: {file_path}")]
    
    sys.path.insert(0, str(Path(__file__).parent.parent.parent / "core"))
    from openqc.parsers import parse_file
    from collections import Counter
    
    result = parse_file(file_path)
    
    if not result.structure:
        return [TextContent(type="text", text="No structure found")]
    
    atoms = result.structure.atoms
    
    # Calculate molecular formula
    element_counts = Counter(atom.element for atom in atoms)
    formula = "".join(
        f"{elem}{count if count > 1 else ''}"
        for elem, count in sorted(element_counts.items())
    )
    
    # Calculate molecular mass (approximate)
    atomic_masses = {
        'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
        'F': 18.998, 'P': 30.974, 'S': 32.065, 'Cl': 35.453,
        'Br': 79.904, 'I': 126.904, 'Fe': 55.845, 'Cu': 63.546,
        'Zn': 65.38, 'Na': 22.990, 'K': 39.098, 'Ca': 40.078,
    }
    
    total_mass = sum(
        atomic_masses.get(atom.element, 0) * count
        for atom.element, count in element_counts.items()
    )
    
    info = {
        "formula": formula,
        "num_atoms": len(atoms),
        "molecular_mass": round(total_mass, 3),
        "elements": dict(element_counts),
        "charge": result.structure.charge,
        "multiplicity": result.structure.multiplicity,
    }
    
    return [TextContent(type="text", text=json.dumps(info, indent=2))]


async def handle_generate_input_template(args: Dict[str, Any]) -> List[TextContent]:
    """Generate input template."""
    software = args["software"].lower()
    job_type = args.get("job_type", "sp")
    functional = args.get("functional", "B3LYP")
    basis = args.get("basis", "def2-SVP")
    
    templates = {
        "gaussian": f"""%nprocshared=8
%mem=16GB
%chk=calculation.chk
#{job_type} {functional}/{basis}

Title

0 1
<Add your molecule coordinates here>

""",
        "vasp": f"""System: {functional} calculation
520.0
<lattice vectors>
<element names>
<atom counts>
Direct
<atomic positions>

""",
        "orca": f"""! {functional} {basis} {job_type} TightSCF
%pal nprocs 8 end
%maxcore 2000

* xyz 0 1
<Add your molecule coordinates here>
*

""",
        "qe": f"""&CONTROL
    calculation = '{job_type}'
    pseudo_dir = './pseudo'
/
&SYSTEM
    ibrav = 0
    nat = <number of atoms>
    ntyp = <number of element types>
    ecutwfc = 60.0
/
&ELECTRONS
/

ATOMIC_SPECIES
<element definitions>

ATOMIC_POSITIONS angstrom
<atomic positions>

CELL_PARAMETERS angstrom
<lattice vectors>

K_POINTS automatic
4 4 4 0 0 0

"""
    }
    
    template = templates.get(software, "Template not available for this software")
    return [TextContent(type="text", text=template)]


async def main():
    """Main entry point."""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(
            read_stream,
            write_stream,
            server.create_initialization_options()
        )


if __name__ == "__main__":
    asyncio.run(main())
