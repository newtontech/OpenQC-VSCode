"""
Parser for ORCA input files.
"""

import re
from typing import List, Optional
from openqc.parsers.base import (
    BaseParser,
    ParsedInput,
    Structure,
    Atom,
    CalculationParams,
)


class ORCAParser(BaseParser):
    """
    Parser for ORCA input files.
    
    Supports ORCA simple input and block input formats.
    """
    
    # Common ORCA functionals
    FUNCTIONALS = {
        'b3lyp', 'pbe0', 'm06', 'm06-2x', 'm06-l', 'wb97x-d',
        'b2plyp', 'mp2', 'dlpno-ccsd(t)', 'hf'
    }
    
    BASIS_SETS = {
        'sto-3g', 'def2-svp', 'def2-tzvp', 'def2-qzvp',
        'cc-pvdz', 'cc-pvtz', 'ma-def2-svp', 'def2-svp/c'
    }
    
    JOB_TYPES = {
        'sp', 'opt', 'freq', 'neb-ts', 'transitionstate',
        'md', 'nmr', 'epr', 'uvvis', 'casscf'
    }
    
    def parse(self, content: str) -> ParsedInput:
        """Parse ORCA input content."""
        self.reset()
        
        lines = content.split('\n')
        
        # Parse components
        simple_input = []
        blocks = {}
        atoms: List[Atom] = []
        charge = 0
        multiplicity = 1
        
        in_block = None
        block_content = []
        in_coords = False
        coord_type = None
        
        for i, line in enumerate(lines):
            stripped = line.strip()
            
            # Skip empty lines and comments
            if not stripped or stripped.startswith('#'):
                continue
            
            # Simple input line (starts with !)
            if stripped.startswith('!'):
                simple_input.extend(stripped[1:].split())
                continue
            
            # Block start
            if stripped.startswith('%'):
                block_match = re.match(r'%(\w+)', stripped)
                if block_match:
                    in_block = block_match.group(1).lower()
                    block_content = [stripped]
                continue
            
            # Block end
            if stripped.lower() == 'end' and in_block:
                block_content.append(stripped)
                blocks[in_block] = '\n'.join(block_content)
                in_block = None
                block_content = []
                continue
            
            # Inside block
            if in_block:
                block_content.append(stripped)
                continue
            
            # Coordinate section start
            if stripped.startswith('*'):
                coord_match = re.match(r'\*\s*(\w+)\s+(\d+)\s+(\d+)', stripped)
                if coord_match:
                    coord_type = coord_match.group(1).lower()
                    charge = int(coord_match.group(2))
                    multiplicity = int(coord_match.group(3))
                    in_coords = True
                continue
            
            # Coordinate section end
            if stripped == '*' and in_coords:
                in_coords = False
                continue
            
            # Parse coordinates
            if in_coords:
                parts = stripped.split()
                if len(parts) >= 4:
                    element = parts[0]
                    try:
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        atoms.append(Atom(
                            element=element,
                            x=x, y=y, z=z,
                            index=len(atoms) + 1
                        ))
                    except ValueError:
                        self.add_error(f"Cannot parse coordinate: {stripped}")
        
        # Parse simple input for parameters
        params = self._parse_simple_input(simple_input)
        params.charge = charge
        params.multiplicity = multiplicity
        
        # Parse blocks for additional parameters
        self._parse_blocks(blocks, params)
        
        # Create structure
        structure = Structure(
            atoms=atoms,
            charge=charge,
            multiplicity=multiplicity
        )
        
        return ParsedInput(
            structure=structure,
            parameters=params,
            raw_content=content,
            file_type="orca",
            errors=self.errors,
            warnings=self.warnings
        )
    
    def _parse_simple_input(self, keywords: List[str]) -> CalculationParams:
        """Parse simple input line keywords."""
        params = CalculationParams()
        
        keywords_lower = [k.lower() for k in keywords]
        
        # Find functional
        for kw in keywords_lower:
            if kw in self.FUNCTIONALS:
                params.functional = kw.upper()
                params.method = kw.upper()
                break
        
        # Find basis set
        for kw in keywords_lower:
            if kw in self.BASIS_SETS or 'def2' in kw or 'cc-p' in kw:
                params.basis = kw.upper()
                break
        
        # Find job type
        for kw in keywords_lower:
            if kw in self.JOB_TYPES:
                params.job_type = kw.upper()
                break
        
        # Check for other keywords
        extras = {}
        if 'tightscf' in keywords_lower:
            extras['scf_convergence'] = 'tight'
        if 'd3bj' in keywords_lower or 'd3zero' in keywords_lower:
            extras['dispersion'] = 'D3'
        if 'cpcm' in keywords_lower or 'smd' in keywords_lower:
            extras['solvation'] = True
        
        params.additional_keywords = extras
        
        return params
    
    def _parse_blocks(self, blocks: dict, params: CalculationParams) -> None:
        """Parse block input for additional parameters."""
        
        # Parse %pal block for parallel settings
        if 'pal' in blocks:
            nprocs_match = re.search(r'nprocs\s+(\d+)', blocks['pal'])
            if nprocs_match:
                params.additional_keywords['nprocs'] = int(nprocs_match.group(1))
        
        # Parse %maxcore
        if 'maxcore' in blocks:
            mem_match = re.search(r'maxcore\s+(\d+)', blocks['maxcore'])
            if mem_match:
                params.additional_keywords['maxcore'] = int(mem_match.group(1))
        
        # Parse %scf block
        if 'scf' in blocks:
            if 'convergence' in blocks['scf'].lower():
                params.additional_keywords['scf_settings'] = 'custom'
        
        # Parse %geom block
        if 'geom' in blocks:
            if 'constraints' in blocks['geom'].lower():
                params.additional_keywords['constraints'] = True
