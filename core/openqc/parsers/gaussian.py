"""
Parser for Gaussian input files (.gjf, .com).
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


class GaussianParser(BaseParser):
    """
    Parser for Gaussian input files.
    
    Supports standard Gaussian input format with:
    - Route section (Link 0 commands)
    - Title section
    - Charge and multiplicity
    - Molecular geometry
    - Additional sections
    """
    
    # Common Gaussian keywords
    FUNCTIONALS = {
        'b3lyp', 'pbe0', 'm06', 'm06-2x', 'm06-l', 'wb97x-d', 
        'cam-b3lyp', 'b2plyp', 'mp2', 'hf', 'ccsd', 'ccsd(t)'
    }
    
    BASIS_SETS = {
        'sto-3g', '3-21g', '6-31g', '6-31g(d)', '6-31g(d,p)',
        '6-311g', '6-311+g(d,p)', 'def2-svp', 'def2-tzvp', 
        'def2-qzvp', 'cc-pvdz', 'cc-pvtz', 'cc-pvqz', 'aug-cc-pvdz'
    }
    
    JOB_TYPES = {
        'sp', 'opt', 'freq', 'irc', 'scan', 'td', 'nmr',
        'force', 'pop', 'density', 'guess', 'stable'
    }
    
    def parse(self, content: str) -> ParsedInput:
        """Parse Gaussian input content."""
        self.reset()
        
        lines = content.split('\n')
        
        # Initialize sections
        route_section = ""
        title = ""
        charge = 0
        multiplicity = 1
        atoms: List[Atom] = []
        
        # Parse sections
        section_idx = 0
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Skip empty lines
            if not line:
                i += 1
                continue
            
            # Route section (starts with % or #)
            if line.startswith('%') or line.startswith('#'):
                route_section = self._parse_route_section(lines, i)
                # Find end of route section
                while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('!'):
                    i += 1
                continue
            
            # Title line
            if section_idx == 0 and not route_section:
                # First non-route line might be title
                if not line[0].isdigit() and not re.match(r'^[A-Za-z]+\s+[\d.-]', line):
                    title = line
                    i += 1
                    section_idx += 1
                    continue
            
            # Charge and multiplicity
            cm_match = re.match(r'^([\d-]+)\s+(\d+)\s*$', line)
            if cm_match:
                charge = int(cm_match.group(1))
                multiplicity = int(cm_match.group(2))
                i += 1
                continue
            
            # Atom line: Element X Y Z or Element Index X Y Z
            atom_match = re.match(
                r'^([A-Za-z]+)(?:\s+\d+)?\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)\s*$',
                line
            )
            if atom_match:
                element = atom_match.group(1).capitalize()
                x = float(atom_match.group(2))
                y = float(atom_match.group(3))
                z = float(atom_match.group(4))
                atoms.append(Atom(element=element, x=x, y=y, z=z, index=len(atoms) + 1))
            
            i += 1
        
        # Parse parameters from route section
        parameters = self._parse_route_parameters(route_section)
        
        # Create structure
        structure = Structure(
            atoms=atoms,
            title=title,
            charge=charge,
            multiplicity=multiplicity
        )
        
        return ParsedInput(
            structure=structure,
            parameters=parameters,
            raw_content=content,
            file_type="gaussian",
            errors=self.errors,
            warnings=self.warnings
        )
    
    def _parse_route_section(self, lines: List[str], start: int) -> str:
        """Parse the route section starting from a given line."""
        route_lines = []
        
        for i in range(start, len(lines)):
            line = lines[i].strip()
            
            if not line or line.startswith('!'):
                break
            
            route_lines.append(line)
        
        return ' '.join(route_lines)
    
    def _parse_route_parameters(self, route: str) -> CalculationParams:
        """Parse calculation parameters from the route section."""
        params = CalculationParams()
        
        # Normalize route
        route_lower = route.lower()
        
        # Find functional
        for func in self.FUNCTIONALS:
            if func in route_lower:
                params.functional = func.upper()
                break
        
        # Find basis set
        for basis in self.BASIS_SETS:
            if basis in route_lower:
                params.basis = basis.upper()
                break
        
        # Find job type
        for job in self.JOB_TYPES:
            if job in route_lower:
                params.job_type = job.upper()
                break
        
        # Method = functional or basis
        if params.functional:
            params.method = params.functional
        
        return params
