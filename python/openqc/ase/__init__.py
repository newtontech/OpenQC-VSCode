"""
OpenQC ASE Integration Module

This module provides ASE (Atomic Simulation Environment) integration
for OpenQC-VSCode, enabling cross-code workflow migration and structure
conversion between different quantum chemistry packages.

Author: NewtonTech
Version: 1.0.0
"""

__version__ = '1.0.0'

from .converter import ASEConverter, ConversionResult
from .utils import validate_atoms, get_supported_formats

__all__ = [
    'ASEConverter',
    'ConversionResult',
    'validate_atoms',
    'get_supported_formats',
]
