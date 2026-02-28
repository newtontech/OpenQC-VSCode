"""
Chemistry utility functions.
"""

import math
from typing import Tuple, List, Dict

# Atomic data
ATOMIC_MASSES: Dict[str, float] = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
    'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
    'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.64, 'As': 74.922, 'Se': 78.96, 'Br': 79.904,
    'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
    'Nb': 92.906, 'Mo': 95.96, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.906,
    'Pd': 106.42, 'Ag': 107.868, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71,
    'Sb': 121.76, 'Te': 127.6, 'I': 126.904, 'Xe': 131.293, 'Cs': 132.905,
    'Ba': 137.327, 'La': 138.905, 'Ce': 140.116, 'Pr': 140.908, 'Nd': 144.242,
    'Pm': 145.0, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.925,
    'Dy': 162.5, 'Ho': 164.930, 'Er': 167.259, 'Tm': 168.934, 'Yb': 173.054,
    'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.948, 'W': 183.84, 'Re': 186.207,
    'Os': 190.23, 'Ir': 192.217, 'Pt': 195.084, 'Au': 196.967, 'Hg': 200.59,
    'Tl': 204.383, 'Pb': 207.2, 'Bi': 208.980, 'Po': 209.0, 'At': 210.0,
    'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.038,
    'Pa': 231.036, 'U': 238.029, 'Np': 237.0, 'Pu': 244.0, 'Am': 243.0,
}

ATOMIC_NUMBERS: Dict[str, int] = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5,
    'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
    'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
    'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
    'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35,
    'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40,
    'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45,
    'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
    'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
}

SYMBOLS_BY_NUMBER = {v: k for k, v in ATOMIC_NUMBERS.items()}


def get_atomic_mass(element: str) -> float:
    """Get atomic mass for an element."""
    element = element.capitalize()
    return ATOMIC_MASSES.get(element, 0.0)


def get_atomic_number(element: str) -> int:
    """Get atomic number for an element."""
    element = element.capitalize()
    return ATOMIC_NUMBERS.get(element, 0)


def get_element_symbol(atomic_number: int) -> str:
    """Get element symbol from atomic number."""
    return SYMBOLS_BY_NUMBER.get(atomic_number, 'X')


def calculate_distance(
    coord1: Tuple[float, float, float],
    coord2: Tuple[float, float, float]
) -> float:
    """Calculate distance between two points in 3D space."""
    return math.sqrt(
        (coord2[0] - coord1[0])**2 +
        (coord2[1] - coord1[1])**2 +
        (coord2[2] - coord1[2])**2
    )


def calculate_angle(
    coord1: Tuple[float, float, float],
    coord2: Tuple[float, float, float],
    coord3: Tuple[float, float, float]
) -> float:
    """Calculate angle in degrees between three points (coord2 is center)."""
    # Vector from coord2 to coord1
    v1 = (
        coord1[0] - coord2[0],
        coord1[1] - coord2[1],
        coord1[2] - coord2[2]
    )
    
    # Vector from coord2 to coord3
    v2 = (
        coord3[0] - coord2[0],
        coord3[1] - coord2[1],
        coord3[2] - coord2[2]
    )
    
    # Calculate dot product and magnitudes
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    
    if mag1 == 0 or mag2 == 0:
        return 0.0
    
    # Calculate angle
    cos_angle = dot / (mag1 * mag2)
    cos_angle = max(-1, min(1, cos_angle))  # Clamp to [-1, 1]
    
    return math.degrees(math.acos(cos_angle))


def calculate_dihedral(
    coord1: Tuple[float, float, float],
    coord2: Tuple[float, float, float],
    coord3: Tuple[float, float, float],
    coord4: Tuple[float, float, float]
) -> float:
    """Calculate dihedral angle in degrees for four points."""
    # Vectors along the bonds
    b1 = (
        coord1[0] - coord2[0],
        coord1[1] - coord2[1],
        coord1[2] - coord2[2]
    )
    b2 = (
        coord2[0] - coord3[0],
        coord2[1] - coord3[1],
        coord2[2] - coord3[2]
    )
    b3 = (
        coord3[0] - coord4[0],
        coord3[1] - coord4[1],
        coord3[2] - coord4[2]
    )
    
    # Normal vectors
    n1 = cross_product(b1, b2)
    n2 = cross_product(b2, b3)
    
    # Normalize
    n1_mag = math.sqrt(n1[0]**2 + n1[1]**2 + n1[2]**2)
    n2_mag = math.sqrt(n2[0]**2 + n2[1]**2 + n2[2]**2)
    
    if n1_mag == 0 or n2_mag == 0:
        return 0.0
    
    n1 = (n1[0]/n1_mag, n1[1]/n1_mag, n1[2]/n1_mag)
    n2 = (n2[0]/n2_mag, n2[1]/n2_mag, n2[2]/n2_mag)
    
    # Calculate angle
    dot = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]
    dot = max(-1, min(1, dot))
    
    return math.degrees(math.acos(dot))


def cross_product(
    v1: Tuple[float, float, float],
    v2: Tuple[float, float, float]
) -> Tuple[float, float, float]:
    """Calculate cross product of two vectors."""
    return (
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    )


def get_center_of_mass(
    coords: List[Tuple[float, float, float]],
    elements: List[str]
) -> Tuple[float, float, float]:
    """Calculate center of mass for a set of atoms."""
    total_mass = 0.0
    com = [0.0, 0.0, 0.0]
    
    for coord, element in zip(coords, elements):
        mass = get_atomic_mass(element)
        total_mass += mass
        com[0] += coord[0] * mass
        com[1] += coord[1] * mass
        com[2] += coord[2] * mass
    
    if total_mass == 0:
        return (0.0, 0.0, 0.0)
    
    return (com[0]/total_mass, com[1]/total_mass, com[2]/total_mass)
