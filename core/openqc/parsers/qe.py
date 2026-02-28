"""
Parser for Quantum ESPRESSO input files.
"""

import re
from typing import Dict, List, Optional, Any
from openqc.parsers.base import (
    BaseParser,
    ParsedInput,
    Structure,
    Atom,
    CalculationParams,
)


class QEParser(BaseParser):
    """
    Parser for Quantum ESPRESSO input files.
    
    Supports pw.x input format with namelists and cards.
    """
    
    def parse(self, content: str) -> ParsedInput:
        """Parse Quantum ESPRESSO input content."""
        self.reset()
        
        # Parse namelists
        namelists = self._parse_namelists(content)
        
        # Parse cards
        cards = self._parse_cards(content)
        
        # Extract structure from ATOMIC_POSITIONS
        atoms: List[Atom] = []
        cell: List[List[float]] = []
        
        if 'ATOMIC_POSITIONS' in cards:
            atoms = self._parse_atomic_positions(cards['ATOMIC_POSITIONS'])
        
        if 'CELL_PARAMETERS' in cards:
            cell = self._parse_cell_parameters(cards['CELL_PARAMETERS'])
        elif 'SYSTEM' in namelists:
            # Use celldm or A, B, C from &SYSTEM
            cell = self._cell_from_system(namelists['SYSTEM'])
        
        # Extract calculation parameters
        params = self._extract_parameters(namelists)
        
        # Create structure
        structure = Structure(
            atoms=atoms,
            cell=cell if cell else None,
            charge=0,
            multiplicity=1
        )
        
        return ParsedInput(
            structure=structure,
            parameters=params,
            raw_content=content,
            file_type="quantumespresso",
            errors=self.errors,
            warnings=self.warnings
        )
    
    def _parse_namelists(self, content: str) -> Dict[str, Dict[str, Any]]:
        """Parse all namelists from the content."""
        namelists = {}
        
        # Find all namelists
        pattern = r'&(\w+)\s*\n(.*?)\s*/'
        matches = re.finditer(pattern, content, re.DOTALL | re.IGNORECASE)
        
        for match in matches:
            name = match.group(1).upper()
            body = match.group(2)
            
            # Parse key-value pairs
            params = {}
            for line in body.split('\n'):
                line = line.strip()
                if not line or line.startswith('!'):
                    continue
                
                # Parse key = value
                kv_match = re.match(r'(\w+)\s*=\s*(.+?)(?:\s*[,\n]|$)', line)
                if kv_match:
                    key = kv_match.group(1).lower()
                    value = kv_match.group(2).strip().rstrip(',')
                    
                    # Convert value type
                    value = self._convert_value(value)
                    params[key] = value
            
            namelists[name] = params
        
        return namelists
    
    def _parse_cards(self, content: str) -> Dict[str, str]:
        """Parse all cards from the content."""
        cards = {}
        
        card_names = [
            'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 
            'CELL_PARAMETERS', 'K_POINTS', 
            'CONSTRAINTS', 'OCCUPATIONS'
        ]
        
        for card_name in card_names:
            pattern = rf'{card_name}\s*(?:\([^)]*\))?\s*\n(.*?)(?={ "|".join(card_names) }|\Z)'
            match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
            if match:
                cards[card_name] = match.group(1).strip()
        
        return cards
    
    def _parse_atomic_positions(self, content: str) -> List[Atom]:
        """Parse ATOMIC_POSITIONS card."""
        atoms = []
        
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            
            parts = line.split()
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
                    self.add_error(f"Cannot parse atomic position: {line}")
        
        return atoms
    
    def _parse_cell_parameters(self, content: str) -> List[List[float]]:
        """Parse CELL_PARAMETERS card."""
        cell = []
        
        for line in content.split('\n'):
            line = line.strip()
            if not line or line.startswith('!'):
                continue
            
            parts = line.split()
            if len(parts) >= 3:
                try:
                    vec = [float(parts[0]), float(parts[1]), float(parts[2])]
                    cell.append(vec)
                except ValueError:
                    self.add_error(f"Cannot parse cell parameter: {line}")
        
        return cell
    
    def _cell_from_system(self, system: Dict[str, Any]) -> List[List[float]]:
        """Construct cell from &SYSTEM parameters."""
        # This is a simplified implementation
        # Full implementation would use celldm or A, B, C, cosAB, etc.
        cell = []
        
        if 'a' in system:
            a = system['a']
            cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        elif 'celldm' in system and isinstance(system['celldm'], list):
            celldm = system['celldm']
            a = celldm[0] if len(celldm) > 0 else 1.0
            cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]
        
        return cell
    
    def _convert_value(self, value: str) -> Any:
        """Convert string value to appropriate type."""
        value = value.strip()
        
        # Boolean
        if value.lower() in ['.true.', 'true']:
            return True
        if value.lower() in ['.false.', 'false']:
            return False
        
        # String (remove quotes)
        if value.startswith('"') and value.endswith('"'):
            return value[1:-1]
        if value.startswith("'") and value.endswith("'"):
            return value[1:-1]
        
        # Number
        try:
            if '.' in value or 'e' in value.lower():
                return float(value)
            return int(value)
        except ValueError:
            pass
        
        return value
    
    def _extract_parameters(self, namelists: Dict[str, Dict]) -> CalculationParams:
        """Extract calculation parameters from namelists."""
        params = CalculationParams()
        
        control = namelists.get('CONTROL', {})
        system = namelists.get('SYSTEM', {})
        electrons = namelists.get('ELECTRONS', {})
        
        # Extract relevant parameters
        params.additional_keywords = {
            'calculation': control.get('calculation', 'scf'),
            'ecutwfc': system.get('ecutwfc'),
            'ecutrho': system.get('ecutrho'),
            'ibrav': system.get('ibrav', 0),
            'nat': system.get('nat'),
            'ntyp': system.get('ntyp'),
            'conv_thr': electrons.get('conv_thr'),
        }
        
        # Job type
        calc = control.get('calculation', 'scf')
        if calc == 'relax':
            params.job_type = 'OPT'
        elif calc == 'md':
            params.job_type = 'MD'
        elif calc == 'vc-relax':
            params.job_type = 'CELL_OPT'
        else:
            params.job_type = 'SP'
        
        return params
