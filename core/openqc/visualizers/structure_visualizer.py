"""
Structure visualization using various backends.
"""

import json
from pathlib import Path
from typing import Optional, Union, Dict, Any, List
from openqc.parsers.base import Structure, Atom


class StructureVisualizer:
    """
    Visualizer for molecular and crystal structures.
    
    Supports multiple backends:
    - 3Dmol.js (web-based, WebGL)
    - PyMOL
    - ASE viewer
    - NGL viewer
    """
    
    def __init__(self, backend: str = '3dmol'):
        """
        Initialize visualizer.
        
        Args:
            backend: Visualization backend ('3dmol', 'pymol', 'ase', 'ngl')
        """
        self.backend = backend.lower()
    
    def visualize(
        self,
        structure: Structure,
        style: str = 'stick',
        **kwargs
    ) -> Any:
        """
        Visualize a structure.
        
        Args:
            structure: Structure to visualize
            style: Visualization style ('stick', 'sphere', 'cartoon', 'line')
            **kwargs: Additional style parameters
            
        Returns:
            Visualization output (depends on backend)
        """
        if self.backend == '3dmol':
            return self._visualize_3dmol(structure, style, **kwargs)
        elif self.backend == 'pymol':
            return self._visualize_pymol(structure, style, **kwargs)
        elif self.backend == 'ase':
            return self._visualize_ase(structure, style, **kwargs)
        else:
            raise ValueError(f"Unknown backend: {self.backend}")
    
    def _visualize_3dmol(
        self, 
        structure: Structure, 
        style: str,
        **kwargs
    ) -> str:
        """Generate 3Dmol.js HTML visualization."""
        # Generate XYZ format for 3Dmol
        xyz_lines = [f"{structure.num_atoms}"]
        xyz_lines.append(structure.title or "OpenQC Structure")
        
        for atom in structure.atoms:
            xyz_lines.append(f"{atom.element} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}")
        
        xyz_data = "\n".join(xyz_lines)
        
        # Generate HTML
        html = self._generate_3dmol_html(xyz_data, style, **kwargs)
        return html
    
    def _generate_3dmol_html(
        self,
        xyz_data: str,
        style: str,
        width: int = 600,
        height: int = 400,
        background: str = 'black',
    ) -> str:
        """Generate complete 3Dmol.js HTML page."""
        return f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.1/3Dmol.min.js"></script>
    <style>
        body {{ margin: 0; padding: 0; }}
        #viewer {{ width: {width}px; height: {height}px; }}
    </style>
</head>
<body>
    <div id="viewer"></div>
    <script>
        var viewer = $3Dmol.createViewer("viewer", {{
            defaultcolors: $3Dmol.rasmolElementColors
        }});
        viewer.setBackgroundColor("{background}");
        
        var xyzData = `{xyz_data}`;
        viewer.addModel(xyzData, "xyz");
        viewer.setStyle({{}}, {style}: {{}});
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>'''
    
    def _visualize_pymol(
        self, 
        structure: Structure,
        style: str,
        **kwargs
    ) -> None:
        """Visualize using PyMOL."""
        try:
            from pymol import cmd
        except ImportError:
            raise ImportError(
                "PyMOL is not installed. Install with:\n"
                "  pip install pymol-open-source\n"
                "  or\n"
                "  conda install -c conda-forge pymol-open-source"
            )
        
        # Create temporary XYZ file
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            f.write(f"{structure.num_atoms}\n")
            f.write(f"{structure.title or 'OpenQC Structure'}\n")
            for atom in structure.atoms:
                f.write(f"{atom.element} {atom.x:.6f} {atom.y:.6f} {atom.z:.6f}\n")
            temp_path = f.name
        
        # Load into PyMOL
        cmd.load(temp_path, "structure")
        
        # Apply style
        if style == 'stick':
            cmd.show('sticks', 'all')
        elif style == 'sphere':
            cmd.show('spheres', 'all')
        elif style == 'cartoon':
            cmd.show('cartoon', 'all')
        elif style == 'line':
            cmd.show('lines', 'all')
    
    def _visualize_ase(
        self, 
        structure: Structure,
        style: str,
        **kwargs
    ) -> None:
        """Visualize using ASE viewer."""
        try:
            from ase import Atoms
            from ase.visualize import view
        except ImportError:
            raise ImportError(
                "ASE is not installed. Install with:\n"
                "  pip install ase"
            )
        
        # Convert to ASE Atoms
        positions = [[atom.x, atom.y, atom.z] for atom in structure.atoms]
        symbols = [atom.element for atom in structure.atoms]
        
        cell = structure.cell if structure.cell else None
        atoms = Atoms(
            symbols=symbols,
            positions=positions,
            cell=cell,
            pbc=cell is not None
        )
        
        view(atoms)
    
    def to_json(self, structure: Structure) -> Dict[str, Any]:
        """Convert structure to JSON format for web visualization."""
        return {
            'atoms': [
                {
                    'element': atom.element,
                    'x': atom.x,
                    'y': atom.y,
                    'z': atom.z,
                    'index': atom.index
                }
                for atom in structure.atoms
            ],
            'cell': structure.cell,
            'title': structure.title,
            'charge': structure.charge,
            'multiplicity': structure.multiplicity,
            'numAtoms': structure.num_atoms
        }


def export_visualization(
    structure: Structure,
    output_file: Union[str, Path],
    format: str = 'png',
    width: int = 1200,
    height: int = 800,
    **kwargs
) -> Path:
    """
    Export structure visualization to file.
    
    Args:
        structure: Structure to visualize
        output_file: Output file path
        format: Output format ('png', 'svg', 'html', 'json')
        width: Image width
        height: Image height
        **kwargs: Additional visualization parameters
        
    Returns:
        Path to the output file
    """
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    if format == 'html':
        viz = StructureVisualizer('3dmol')
        html = viz._visualize_3dmol(structure, **kwargs)
        output_path.write_text(html)
    
    elif format == 'json':
        viz = StructureVisualizer('3dmol')
        data = viz.to_json(structure)
        output_path.write_text(json.dumps(data, indent=2))
    
    elif format in ['png', 'svg']:
        # Use PyMOL for image export if available
        try:
            viz = StructureVisualizer('pymol')
            viz._visualize_pymol(structure, **kwargs)
            
            from pymol import cmd
            if format == 'png':
                cmd.png(str(output_path), width=width, height=height)
            else:
                cmd.svg(str(output_path))
        except ImportError:
            # Fall back to generating HTML
            viz = StructureVisualizer('3dmol')
            html = viz._visualize_3dmol(structure, **kwargs)
            html_path = output_path.with_suffix('.html')
            html_path.write_text(html)
            raise RuntimeError(
                f"Cannot export to {format} without PyMOL. "
                f"HTML visualization saved to {html_path}"
            )
    
    return output_path
