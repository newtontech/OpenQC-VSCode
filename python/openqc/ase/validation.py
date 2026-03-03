"""
ASE Validation Utilities

Backend validation utilities for structure comparison,
energy/force validation, and round-trip conversion testing.
"""

from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import json
import sys
from dataclasses import dataclass, asdict
from ase import Atoms
from ase.io import read, write
import numpy as np


@dataclass
class ValidationResult:
    """Result of a validation operation."""
    success: bool
    test_name: str
    errors: List[str]
    warnings: List[str]
    metrics: Dict[str, Any]
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class StructureComparison:
    """Result of structure comparison."""
    positions_match: bool
    cell_match: bool
    pbc_match: bool
    elements_match: bool
    natoms_match: bool
    max_position_diff: float
    max_cell_diff: float
    rms_position_diff: float
    rms_cell_diff: float
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class EnergyForceComparison:
    """Result of energy/force comparison."""
    energy_match: bool
    forces_match: bool
    energy_diff: float
    max_force_diff: float
    rms_force_diff: float
    relative_energy_diff: float
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class MigrationValidator:
    """
    Validation utilities for migration between quantum chemistry codes.
    """
    
    def __init__(self, tolerance: float = 1e-6, energy_tolerance: float = 1e-3):
        """
        Initialize validator.
        
        Parameters
        ----------
        tolerance : float
            Position/cell tolerance in Angstrom
        energy_tolerance : float
            Energy tolerance in eV
        """
        self.tolerance = tolerance
        self.energy_tolerance = energy_tolerance
    
    def compare_structures(
        self,
        atoms1: Atoms,
        atoms2: Atoms,
        ignore_order: bool = True
    ) -> StructureComparison:
        """
        Compare two ASE Atoms structures.
        
        Parameters
        ----------
        atoms1 : Atoms
            First structure
        atoms2 : Atoms
            Second structure
        ignore_order : bool
            Whether to ignore atom ordering
        
        Returns
        -------
        StructureComparison
            Comparison result
        """
        # Check atom count
        natoms_match = len(atoms1) == len(atoms2)
        
        # Check elements
        symbols1 = atoms1.get_chemical_symbols()
        symbols2 = atoms2.get_chemical_symbols()
        
        if ignore_order:
            elements_match = sorted(symbols1) == sorted(symbols2)
        else:
            elements_match = symbols1 == symbols2
        
        # Check PBC
        pbc_match = np.all(atoms1.pbc == atoms2.pbc)
        
        # Check cell
        cell1 = atoms1.get_cell()
        cell2 = atoms2.get_cell()
        
        if atoms1.pbc.any() and atoms2.pbc.any():
            cell_diff = np.abs(cell1 - cell2)
            max_cell_diff = np.max(cell_diff)
            rms_cell_diff = np.sqrt(np.mean(cell_diff**2))
            cell_match = max_cell_diff < self.tolerance
        else:
            max_cell_diff = 0.0
            rms_cell_diff = 0.0
            cell_match = True
        
        # Check positions
        pos1 = atoms1.get_positions()
        pos2 = atoms2.get_positions()
        
        if ignore_order:
            # Sort by position for comparison
            idx1 = np.lexsort((pos1[:, 2], pos1[:, 1], pos1[:, 0]))
            idx2 = np.lexsort((pos2[:, 2], pos2[:, 1], pos2[:, 0]))
            pos1 = pos1[idx1]
            pos2 = pos2[idx2]
        
        pos_diff = np.abs(pos1 - pos2)
        max_position_diff = np.max(pos_diff)
        rms_position_diff = np.sqrt(np.mean(pos_diff**2))
        positions_match = max_position_diff < self.tolerance
        
        return StructureComparison(
            positions_match=positions_match,
            cell_match=cell_match,
            pbc_match=pbc_match,
            elements_match=elements_match,
            natoms_match=natoms_match,
            max_position_diff=max_position_diff,
            max_cell_diff=max_cell_diff,
            rms_position_diff=rms_position_diff,
            rms_cell_diff=rms_cell_diff
        )
    
    def compare_energy_forces(
        self,
        energy1: float,
        forces1: np.ndarray,
        energy2: float,
        forces2: np.ndarray
    ) -> EnergyForceComparison:
        """
        Compare energy and forces between two calculations.
        
        Parameters
        ----------
        energy1 : float
            First energy in eV
        forces1 : np.ndarray
            First forces in eV/Angstrom
        energy2 : float
            Second energy in eV
        forces2 : np.ndarray
            Second forces in eV/Angstrom
        
        Returns
        -------
        EnergyForceComparison
            Comparison result
        """
        # Energy comparison
        energy_diff = abs(energy1 - energy2)
        relative_energy_diff = energy_diff / max(abs(energy1), abs(energy2), 1e-10)
        energy_match = energy_diff < self.energy_tolerance
        
        # Force comparison
        if forces1.shape == forces2.shape:
            force_diff = np.abs(forces1 - forces2)
            max_force_diff = np.max(force_diff)
            rms_force_diff = np.sqrt(np.mean(force_diff**2))
            forces_match = max_force_diff < self.energy_tolerance * 10  # Forces typically larger
        else:
            max_force_diff = float('inf')
            rms_force_diff = float('inf')
            forces_match = False
        
        return EnergyForceComparison(
            energy_match=energy_match,
            forces_match=forces_match,
            energy_diff=energy_diff,
            max_force_diff=max_force_diff,
            rms_force_diff=rms_force_diff,
            relative_energy_diff=relative_energy_diff
        )
    
    def validate_round_trip(
        self,
        original_file: str,
        intermediate_format: str,
        output_file: str,
        original_format: str = None
    ) -> ValidationResult:
        """
        Validate round-trip conversion.
        
        Parameters
        ----------
        original_file : str
            Path to original file
        intermediate_format : str
            Intermediate ASE format
        output_file : str
            Path to output file
        original_format : str, optional
            Original file format
        
        Returns
        -------
        ValidationResult
            Validation result
        """
        errors = []
        warnings = []
        metrics = {}
        
        try:
            # Read original
            atoms_original = read(original_file, format=original_format)
            if atoms_original is None:
                return ValidationResult(
                    success=False,
                    test_name="round_trip",
                    errors=[f"Failed to read original file: {original_file}"],
                    warnings=[],
                    metrics={}
                )
            
            # Write intermediate
            write(output_file, atoms_original, format=intermediate_format)
            
            # Read back
            atoms_round_trip = read(output_file, format=intermediate_format)
            if atoms_round_trip is None:
                return ValidationResult(
                    success=False,
                    test_name="round_trip",
                    errors=[f"Failed to read round-trip file: {output_file}"],
                    warnings=[],
                    metrics={}
                )
            
            # Compare
            comparison = self.compare_structures(atoms_original, atoms_round_trip)
            
            metrics["natoms"] = len(atoms_original)
            metrics["max_position_diff"] = comparison.max_position_diff
            metrics["rms_position_diff"] = comparison.rms_position_diff
            metrics["max_cell_diff"] = comparison.max_cell_diff
            metrics["rms_cell_diff"] = comparison.rms_cell_diff
            
            if not comparison.natoms_match:
                errors.append(f"Atom count mismatch: {len(atoms_original)} vs {len(atoms_round_trip)}")
            
            if not comparison.elements_match:
                errors.append("Element mismatch")
            
            if not comparison.positions_match:
                errors.append(f"Position mismatch: max diff = {comparison.max_position_diff:.6f} Angstrom")
            
            if not comparison.cell_match:
                errors.append(f"Cell mismatch: max diff = {comparison.max_cell_diff:.6f} Angstrom")
            
            if not comparison.pbc_match:
                warnings.append("PBC mismatch")
            
            success = len(errors) == 0
            
        except Exception as e:
            errors.append(str(e))
            success = False
        
        return ValidationResult(
            success=success,
            test_name="round_trip",
            errors=errors,
            warnings=warnings,
            metrics=metrics
        )
    
    def validate_migration(
        self,
        source_file: str,
        target_file: str,
        source_format: str = None,
        target_format: str = None
    ) -> ValidationResult:
        """
        Validate migration between formats.
        
        Parameters
        ----------
        source_file : str
            Path to source file
        target_file : str
            Path to target file
        source_format : str, optional
            Source format
        target_format : str, optional
            Target format
        
        Returns
        -------
        ValidationResult
            Validation result
        """
        errors = []
        warnings = []
        metrics = {}
        
        try:
            # Read source
            atoms_source = read(source_file, format=source_format)
            if atoms_source is None:
                return ValidationResult(
                    success=False,
                    test_name="migration",
                    errors=[f"Failed to read source file: {source_file}"],
                    warnings=[],
                    metrics={}
                )
            
            # Read target
            atoms_target = read(target_file, format=target_format)
            if atoms_target is None:
                return ValidationResult(
                    success=False,
                    test_name="migration",
                    errors=[f"Failed to read target file: {target_file}"],
                    warnings=[],
                    metrics={}
                )
            
            # Compare
            comparison = self.compare_structures(atoms_source, atoms_target)
            
            metrics["source_natoms"] = len(atoms_source)
            metrics["target_natoms"] = len(atoms_target)
            metrics["max_position_diff"] = comparison.max_position_diff
            metrics["rms_position_diff"] = comparison.rms_position_diff
            metrics["max_cell_diff"] = comparison.max_cell_diff
            
            if not comparison.natoms_match:
                errors.append(f"Atom count mismatch: {len(atoms_source)} vs {len(atoms_target)}")
            
            if not comparison.elements_match:
                errors.append("Element mismatch")
            
            if not comparison.positions_match:
                warnings.append(f"Position difference: max = {comparison.max_position_diff:.6f} Angstrom")
            
            if not comparison.cell_match:
                warnings.append(f"Cell difference: max = {comparison.max_cell_diff:.6f} Angstrom")
            
            success = len(errors) == 0
            
        except Exception as e:
            errors.append(str(e))
            success = False
        
        return ValidationResult(
            success=success,
            test_name="migration",
            errors=errors,
            warnings=warnings,
            metrics=metrics
        )
    
    def generate_validation_report(
        self,
        results: List[ValidationResult]
    ) -> Dict[str, Any]:
        """
        Generate a validation report from multiple results.
        
        Parameters
        ----------
        results : List[ValidationResult]
            List of validation results
        
        Returns
        -------
        Dict[str, Any]
            Validation report
        """
        total_tests = len(results)
        passed_tests = sum(1 for r in results if r.success)
        failed_tests = total_tests - passed_tests
        
        all_errors = []
        all_warnings = []
        
        for r in results:
            all_errors.extend(r.errors)
            all_warnings.extend(r.warnings)
        
        return {
            "summary": {
                "total_tests": total_tests,
                "passed": passed_tests,
                "failed": failed_tests,
                "pass_rate": passed_tests / total_tests if total_tests > 0 else 0
            },
            "errors": all_errors,
            "warnings": all_warnings,
            "details": [r.to_dict() for r in results]
        }


def main():
    """CLI interface for validation."""
    if len(sys.argv) < 2:
        print("Usage: python validation.py <command> [args]")
        print("Commands:")
        print("  compare <file1> <file2>")
        print("  roundtrip <original> <intermediate_format> <output>")
        print("  migrate <source> <target>")
        sys.exit(1)
    
    command = sys.argv[1]
    validator = MigrationValidator()
    
    if command == "compare":
        if len(sys.argv) < 4:
            print("Error: two files required")
            sys.exit(1)
        
        atoms1 = read(sys.argv[2])
        atoms2 = read(sys.argv[3])
        
        result = validator.compare_structures(atoms1, atoms2)
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == "roundtrip":
        if len(sys.argv) < 5:
            print("Error: original, format, and output required")
            sys.exit(1)
        
        result = validator.validate_round_trip(
            sys.argv[2],
            sys.argv[3],
            sys.argv[4]
        )
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == "migrate":
        if len(sys.argv) < 4:
            print("Error: source and target files required")
            sys.exit(1)
        
        result = validator.validate_migration(sys.argv[2], sys.argv[3])
        print(json.dumps(result.to_dict(), indent=2))
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)


if __name__ == "__main__":
    main()
