"""
OpenQC Command Line Interface
"""

import click
from pathlib import Path
from rich.console import Console
from rich.table import Table

console = Console()


@click.group()
@click.version_option(version="0.1.0")
def main():
    """OpenQC - Quantum Chemistry Input File Toolkit"""
    pass


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--format", "-f", help="Output format")
@click.option("--output", "-o", help="Output file path")
def convert(input_file: str, format: str, output: str):
    """Convert between quantum chemistry file formats."""
    from openqc.converters import convert_format
    
    input_path = Path(input_file)
    
    if output is None:
        output = str(input_path.with_suffix(f".{format}" if format else ".xyz"))
    
    try:
        result = convert_format(input_path, output, format)
        console.print(f"[green]✓[/green] Converted to: {result}")
    except Exception as e:
        console.print(f"[red]✗[/red] Conversion failed: {e}")


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
def parse(input_file: str):
    """Parse and display QC input file information."""
    from openqc.parsers import parse_file
    import json
    
    result = parse_file(input_file)
    
    # Display results
    table = Table(title=f"File: {input_file}")
    table.add_column("Property", style="cyan")
    table.add_column("Value", style="green")
    
    table.add_row("File Type", result.file_type)
    
    if result.structure:
        table.add_row("Num Atoms", str(result.structure.num_atoms))
        table.add_row("Elements", ", ".join(set(result.structure.elements)))
        table.add_row("Charge", str(result.structure.charge))
        table.add_row("Multiplicity", str(result.structure.multiplicity))
    
    if result.parameters:
        if result.parameters.functional:
            table.add_row("Functional", result.parameters.functional)
        if result.parameters.basis:
            table.add_row("Basis Set", result.parameters.basis)
        if result.parameters.job_type:
            table.add_row("Job Type", result.parameters.job_type)
    
    console.print(table)
    
    if result.errors:
        console.print("\n[red]Errors:[/red]")
        for error in result.errors:
            console.print(f"  • {error}")
    
    if result.warnings:
        console.print("\n[yellow]Warnings:[/yellow]")
        for warning in result.warnings:
            console.print(f"  • {warning}")


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option("--style", "-s", default="stick", help="Visualization style")
@click.option("--output", "-o", help="Output HTML file")
def visualize(input_file: str, style: str, output: str):
    """Generate 3D visualization of structure."""
    from openqc.parsers import parse_file
    from openqc.visualizers import StructureVisualizer, export_visualization
    
    result = parse_file(input_file)
    
    if not result.structure:
        console.print("[red]✗[/red] No structure found in file")
        return
    
    if output is None:
        output = str(Path(input_file).with_suffix(".html"))
    
    viz = StructureVisualizer("3dmol")
    html = viz.visualize(result.structure, style=style)
    
    Path(output).write_text(html)
    console.print(f"[green]✓[/green] Visualization saved to: {output}")


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
def validate(input_file: str):
    """Validate a QC input file."""
    from openqc.parsers import parse_file
    
    result = parse_file(input_file)
    
    if result.errors:
        console.print("[red]✗ Validation failed[/red]")
        for error in result.errors:
            console.print(f"  • {error}")
    else:
        console.print("[green]✓ Validation passed[/green]")
    
    if result.warnings:
        console.print("\n[yellow]Warnings:[/yellow]")
        for warning in result.warnings:
            console.print(f"  • {warning}")


@main.command()
@click.option("--software", type=click.Choice(["gaussian", "vasp", "orca", "qe"]), required=True)
@click.option("--job-type", default="sp", help="Calculation type")
@click.option("--functional", default="B3LYP", help="DFT functional")
@click.option("--basis", default="def2-SVP", help="Basis set")
@click.option("--output", "-o", help="Output file")
def template(software: str, job_type: str, functional: str, basis: str, output: str):
    """Generate input file template."""
    templates = {
        "gaussian": f"""%nprocshared=8
%mem=16GB
%chk=calculation.chk
#{job_type} {functional}/{basis}

Title

0 1
<Add coordinates here>

""",
        "vasp": f"""System
1.0
<lattice vectors>
<elements>
<atom counts>
Direct
<positions>

""",
        "orca": f"""! {functional} {basis} {job_type} TightSCF
%pal nprocs 8 end
%maxcore 2000

* xyz 0 1
<Add coordinates here>
*

""",
        "qe": f"""&CONTROL
    calculation = '{job_type}'
/
&SYSTEM
    ibrav = 0
    nat = <N>
    ntyp = <N>
    ecutwfc = 60.0
/
&ELECTRONS
/

ATOMIC_POSITIONS angstrom
<positions>

CELL_PARAMETERS angstrom
<lattice>

K_POINTS automatic
4 4 4 0 0 0

"""
    }
    
    content = templates.get(software, "")
    
    if output:
        Path(output).write_text(content)
        console.print(f"[green]✓[/green] Template saved to: {output}")
    else:
        console.print(content)


@main.command()
def formats():
    """List supported file formats."""
    console.print("\n[bold]Supported Input Formats:[/bold]")
    
    formats_table = Table()
    formats_table.add_column("Software", style="cyan")
    formats_table.add_column("Extensions", style="green")
    formats_table.add_column("Status", style="yellow")
    
    formats_table.add_row("Gaussian", ".gjf, .com", "✓ Full")
    formats_table.add_row("VASP", "POSCAR, CONTCAR, INCAR", "✓ Full")
    formats_table.add_row("Quantum ESPRESSO", ".pw, .scf, .nscf", "✓ Full")
    formats_table.add_row("ORCA", ".inp", "✓ Full")
    formats_table.add_row("CP2K", ".inp", "✓ Full")
    formats_table.add_row("XYZ", ".xyz", "✓ Full")
    formats_table.add_row("CIF", ".cif", "✓ Full")
    formats_table.add_row("PDB", ".pdb", "✓ Full")
    
    console.print(formats_table)


if __name__ == "__main__":
    main()
