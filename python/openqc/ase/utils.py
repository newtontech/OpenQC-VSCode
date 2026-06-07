"""
ASE Utility Functions

Helper functions for ASE Atoms validation and format handling.
"""

from typing import List, Dict, Any, Optional
from ase import Atoms
import numpy as np


def validate_atoms(atoms: Atoms) -> Dict[str, Any]:
    """
    Validate ASE Atoms object and return validation results.
    
    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object to validate
    
    Returns
    -------
    Dict[str, Any]
        Validation results with keys:
        - valid: bool - Whether the structure is valid
        - errors: List[str] - List of validation errors
        - warnings: List[str] - List of validation warnings
        - info: Dict[str, Any] - Structure information
    """
    errors = []
    warnings = []
    info = {}
    
    # Check if atoms object is valid
    if not isinstance(atoms, Atoms):
        errors.append("Input is not an ASE Atoms object")
        return {'valid': False, 'errors': errors, 'warnings': warnings, 'info': info}
    
    # Check number of atoms
    natoms = len(atoms)
    if natoms == 0:
        errors.append("Structure has no atoms")
    else:
        info['natoms'] = natoms
    
    # Check atomic positions
    positions = atoms.get_positions()
    if np.any(np.isnan(positions)):
        errors.append("Structure contains NaN positions")
    
    if np.any(np.isinf(positions)):
        errors.append("Structure contains infinite positions")
    
    # Check for very large coordinates (might indicate unit cell issues)
    max_coord = np.max(np.abs(positions))
    if max_coord > 1000:
        warnings.append(f"Very large coordinates detected (max: {max_coord:.2f}). Check unit cell.")
    
    # Check cell (if periodic)
    if atoms.pbc.any():
        cell = atoms.get_cell()
        if np.linalg.det(cell) == 0:
            errors.append("Unit cell has zero volume")
        else:
            info['cell_volume'] = np.linalg.det(cell)
            info['periodic'] = True
    else:
        info['periodic'] = False
    
    # Check atomic numbers
    atomic_numbers = atoms.get_atomic_numbers()
    if np.any(atomic_numbers <= 0) or np.any(atomic_numbers > 118):
        errors.append("Invalid atomic numbers detected")
    
    # Get element information
    symbols = atoms.get_chemical_symbols()
    unique_symbols = list(set(symbols))
    info['elements'] = unique_symbols
    info['formula'] = atoms.get_chemical_formula()
    
    # Check for constraints
    if atoms.constraints:
        info['constraints'] = len(atoms.constraints)
        warnings.append(f"Structure has {len(atoms.constraints)} constraint(s)")
    
    return {
        'valid': len(errors) == 0,
        'errors': errors,
        'warnings': warnings,
        'info': info
    }


def get_supported_formats() -> Dict[str, Dict[str, Any]]:
    """
    Get dictionary of supported formats for ASE conversion.
    
    Returns
    -------
    Dict[str, Dict[str, Any]]
        Dictionary with format names as keys and metadata as values
    """
    return {
        'vasp': {
            'name': 'VASP',
            'extensions': ['POSCAR', 'CONTCAR'],
            'read': True,
            'write': True,
            'periodic': True,
            'description': 'VASP structure format'
        },
        'cp2k': {
            'name': 'CP2K',
            'extensions': ['.inp'],
            'read': True,
            'write': True,
            'periodic': True,
            'description': 'CP2K input format'
        },
        'qe': {
            'name': 'Quantum ESPRESSO',
            'extensions': ['.in', '.pw.in'],
            'read': True,
            'write': True,
            'periodic': True,
            'description': 'Quantum ESPRESSO input format'
        },
        'gaussian': {
            'name': 'Gaussian',
            'extensions': ['.com', '.gjf'],
            'read': True,
            'write': True,
            'periodic': False,
            'description': 'Gaussian input format'
        },
        'orca': {
            'name': 'ORCA',
            'extensions': ['.inp'],
            'read': True,
            'write': True,
            'periodic': False,
            'description': 'ORCA input format'
        },
        'nwchem': {
            'name': 'NWChem',
            'extensions': ['.nw', '.nwinp'],
            'read': True,
            'write': True,
            'periodic': False,
            'description': 'NWChem input format'
        },
        'gamess': {
            'name': 'GAMESS',
            'extensions': ['.inp'],
            'read': True,
            'write': True,
            'periodic': False,
            'description': 'GAMESS input format'
        },
        'lammps': {
            'name': 'LAMMPS',
            'extensions': ['.lmp', '.data'],
            'read': True,
            'write': True,
            'periodic': True,
            'description': 'LAMMPS data format'
        },
        'xyz': {
            'name': 'XYZ',
            'extensions': ['.xyz'],
            'read': True,
            'write': True,
            'periodic': False,
            'description': 'Generic XYZ format'
        },
        'pdb': {
            'name': 'PDB',
            'extensions': ['.pdb'],
            'read': True,
            'write': True,
            'periodic': False,
            'description': 'Protein Data Bank format'
        },
        'cif': {
            'name': 'CIF',
            'extensions': ['.cif'],
            'read': True,
            'write': True,
            'periodic': True,
            'description': 'Crystallographic Information File'
        }
    }


def estimate_kpoints_from_density(
    atoms: Atoms,
    kpt_density: float = 3.0,
    even: bool = True
) -> tuple:
    """
    Estimate k-point grid from density.
    
    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object with unit cell
    kpt_density : float
        K-point density (points per Å⁻¹)
    even : bool
        Whether to make grid even numbers
    
    Returns
    -------
    tuple
        (kx, ky, kz) k-point grid
    """
    if not atoms.pbc.any():
        return (1, 1, 1)
    
    cell = atoms.get_cell()
    rec_cell = 2 * np.pi * np.linalg.inv(cell).T
    kpt_grid = []
    
    for i in range(3):
        if atoms.pbc[i]:
            length = np.linalg.norm(rec_cell[i])
            kpt = max(1, int(length * kpt_density))
            if even and kpt % 2 != 0:
                kpt += 1
            kpt_grid.append(kpt)
        else:
            kpt_grid.append(1)
    
    return tuple(kpt_grid)


def suggest_pseudopotentials(
    atoms: Atoms,
    code: str = 'vasp',
    functional: str = 'PBE'
) -> Dict[str, str]:
    """
    Suggest pseudopotentials for each element.
    
    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object
    code : str
        Target quantum chemistry code
    functional : str
        Exchange-correlation functional
    
    Returns
    -------
    Dict[str, str]
        Dictionary mapping element symbols to pseudopotential names
    """
    symbols = atoms.get_chemical_symbols()
    unique_symbols = list(set(symbols))
    
    suggestions = {}
    
    if code.lower() == 'vasp':
        # VASP POTCAR suggestions
        for symbol in unique_symbols:
            if functional.upper() == 'PBE':
                suggestions[symbol] = f"{symbol}_pv" if symbol in ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'] else symbol
            elif functional.upper() == 'LDA':
                suggestions[symbol] = symbol
    elif code.lower() == 'cp2k':
        # CP2K GTH pseudopotential suggestions
        for symbol in unique_symbols:
            suggestions[symbol] = f"GTH-{functional}"
    elif code.lower() == 'qe':
        # Quantum ESPRESSO suggestions
        for symbol in unique_symbols:
            suggestions[symbol] = f"{symbol}.{functional}.UPF"
    
    return suggestions
