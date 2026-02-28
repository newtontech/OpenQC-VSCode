"""
Parser for VASP input files (POSCAR, CONTCAR, INCAR).
"""

import re
from typing import List, Optional, Tuple
from openqc.parsers.base import (
    BaseParser,
    ParsedInput,
    Structure,
    Atom,
    CalculationParams,
)


class VASPParser(BaseParser):
    """
    Parser for VASP POSCAR/CONTCAR files.
    
    Supports both standard and VASP 5+ formats with element names.
    """
    
    def parse(self, content: str) -> ParsedInput:
        """Parse VASP POSCAR/CONTCAR content."""
        self.reset()
        
        lines = [l for l in content.split('\n')]
        
        if len(lines) < 8:
            self.add_error("File too short to be a valid POSCAR")
            return ParsedInput(
                raw_content=content,
                file_type="vasp",
                errors=self.errors,
                warnings=self.warnings
            )
        
        # Line 0: Comment/title
        title = lines[0].strip()
        
        # Line 1: Scaling factor
        try:
            scaling = float(lines[1].strip())
        except ValueError:
            self.add_error(f"Invalid scaling factor: {lines[1]}")
            scaling = 1.0
        
        # Lines 2-4: Lattice vectors
        lattice = []
        for i in range(2, 5):
            try:
                vec = [float(x) for x in lines[i].split()]
                if len(vec) != 3:
                    self.add_error(f"Invalid lattice vector on line {i+1}")
                lattice.append([v * scaling for v in vec])
            except ValueError:
                self.add_error(f"Cannot parse lattice vector on line {i+1}")
        
        # Line 5: Element names (VASP 5+) or atom counts (VASP 4)
        line5 = lines[5].strip().split()
        
        # Check if this is VASP 5+ format (element names first)
        elements = []
        try:
            # Try to parse as numbers (VASP 4)
            counts = [int(x) for x in line5]
            element_idx = 6
        except ValueError:
            # VASP 5+ format: element names on line 5, counts on line 6
            elements = line5
            try:
                counts = [int(x) for x in lines[6].strip().split()]
            except ValueError:
                self.add_error("Cannot parse atom counts")
                counts = []
            element_idx = 7
        
        # Line 6 or 7: Coordinate type
        coord_line = lines[element_idx].strip().lower() if element_idx < len(lines) else ""
        
        # Determine coordinate type
        selective_dynamics = False
        if coord_line.startswith('s'):
            selective_dynamics = True
            element_idx += 1
            coord_type = lines[element_idx].strip().lower() if element_idx < len(lines) else "direct"
        else:
            coord_type = coord_line
        
        direct_coords = not coord_type.startswith('c') and not coord_type.startswith('k')
        element_idx += 1
        
        # Parse atomic positions
        atoms: List[Atom] = []
        atom_idx = 0
        
        for i, count in enumerate(counts):
            element = elements[i] if i < len(elements) else f"X{i+1}"
            
            for j in range(count):
                if element_idx >= len(lines):
                    self.add_warning(f"Unexpected end of file while reading atoms")
                    break
                
                pos_line = lines[element_idx].strip().split()
                if not pos_line:
                    element_idx += 1
                    continue
                
                try:
                    x, y, z = float(pos_line[0]), float(pos_line[1]), float(pos_line[2])
                    
                    # Convert direct to Cartesian if needed
                    if direct_coords and lattice:
                        cart = self._direct_to_cartesian(x, y, z, lattice)
                        x, y, z = cart
                    
                    atoms.append(Atom(
                        element=element,
                        x=x * scaling if not direct_coords else x,
                        y=y * scaling if not direct_coords else y,
                        z=z * scaling if not direct_coords else z,
                        index=atom_idx + 1
                    ))
                    
                except (ValueError, IndexError):
                    self.add_error(f"Cannot parse atom position on line {element_idx + 1}")
                
                element_idx += 1
                atom_idx += 1
        
        # Create structure
        structure = Structure(
            atoms=atoms,
            cell=lattice if lattice else None,
            title=title,
        )
        
        return ParsedInput(
            structure=structure,
            parameters=CalculationParams(),
            raw_content=content,
            file_type="vasp",
            errors=self.errors,
            warnings=self.warnings
        )
    
    def _direct_to_cartesian(
        self, 
        u: float, 
        v: float, 
        w: float, 
        lattice: List[List[float]]
    ) -> Tuple[float, float, float]:
        """Convert direct (fractional) coordinates to Cartesian."""
        x = u * lattice[0][0] + v * lattice[1][0] + w * lattice[2][0]
        y = u * lattice[0][1] + v * lattice[1][1] + w * lattice[2][1]
        z = u * lattice[0][2] + v * lattice[1][2] + w * lattice[2][2]
        return (x, y, z)
    
    def parse_incar(self, content: str) -> CalculationParams:
        """Parse VASP INCAR file."""
        params = CalculationParams()
        keywords = {}
        
        for line in content.split('\n'):
            line = line.strip()
            
            # Skip comments and empty lines
            if not line or line.startswith('#') or line.startswith('!'):
                continue
            
            # Parse key = value
            match = re.match(r'(\w+)\s*=\s*(.+)', line)
            if match:
                key = match.group(1).upper()
                value = match.group(2).strip()
                
                # Remove comments from value
                if '#' in value:
                    value = value.split('#')[0].strip()
                if '!' in value:
                    value = value.split('!')[0].strip()
                
                # Convert value type
                if value.lower() in ['.true.', 'true', 't']:
                    value = True
                elif value.lower() in ['.false.', 'false', 'f']:
                    value = False
                else:
                    try:
                        value = float(value) if '.' in value else int(value)
                    except ValueError:
                        pass
                
                keywords[key] = value
        
        params.additional_keywords = keywords
        return params
