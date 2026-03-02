/**
 * Structure converter for various quantum chemistry formats
 *
 * Phase 2: 3D Visualization - Format Conversion
 *
 * This module converts parsed quantum chemistry files into
 * the format expected by ThreeJsRenderer.
 */

import { VASPParser, POSCARData } from '../parsers/VASPParser';
import { Atom, MolecularStructure } from './types';
import { Molecule3D } from './Molecule3D';

export interface ParsedStructure {
  atoms: Atom[];
  lattice?: number[][];
  name?: string;
}

/**
 * Convert various quantum chemistry formats to MolecularStructure
 */
export class StructureConverter {
  private molecule3D: Molecule3D;

  constructor() {
    this.molecule3D = new Molecule3D();
  }

  /**
   * Convert VASP POSCAR content to MolecularStructure
   *
   * @param content - POSCAR file content
   * @param filename - Optional filename
   * @returns MolecularStructure for rendering
   */
  public fromPOSCAR(content: string, filename?: string): MolecularStructure {
    const parser = new VASPParser(content, 'POSCAR');
    const result = parser.parseInput();

    if (result.errors.length > 0) {
      throw new Error(`Failed to parse POSCAR: ${result.errors[0].message}`);
    }

    // Extract atoms using parser methods
    const atomTypes = parser.getAtomTypes();
    const atomCounts = parser.getAtomCounts();
    const coordinates = parser.getCoordinates();
    const latticeVectors = parser.getLatticeVectors();

    // Convert to Atom format
    const atoms: Atom[] = [];
    let atomIndex = 0;

    for (let i = 0; i < atomTypes.length; i++) {
      const element = atomTypes[i];
      const count = atomCounts[i];

      for (let j = 0; j < count && atomIndex < coordinates.length; j++) {
        const coord = coordinates[atomIndex];
        atoms.push({
          element,
          x: coord[0],
          y: coord[1],
          z: coord[2],
        });
        atomIndex++;
      }
    }

    // Convert lattice vectors if present
    const lattice =
      latticeVectors.length > 0
        ? {
            a: latticeVectors[0] as [number, number, number],
            b: latticeVectors[1] as [number, number, number],
            c: latticeVectors[2] as [number, number, number],
          }
        : undefined;

    return {
      atoms,
      lattice,
      name: filename || 'POSCAR',
    };
  }

  /**
   * Convert Gaussian input file content to MolecularStructure
   *
   * @param content - Gaussian input file content
   * @param filename - Optional filename
   * @returns MolecularStructure for rendering
   */
  public fromGaussian(content: string, filename?: string): MolecularStructure {
    const atoms = this.molecule3D.parseAtoms(content, 'Gaussian');

    return {
      atoms,
      name: filename || 'Gaussian Input',
    };
  }

  /**
   * Convert ORCA input file content to MolecularStructure
   *
   * @param content - ORCA input file content
   * @param filename - Optional filename
   * @returns MolecularStructure for rendering
   */
  public fromORCA(content: string, filename?: string): MolecularStructure {
    const atoms = this.molecule3D.parseAtoms(content, 'ORCA');

    return {
      atoms,
      name: filename || 'ORCA Input',
    };
  }

  /**
   * Convert CP2K input file content to MolecularStructure
   *
   * @param content - CP2K input file content
   * @param filename - Optional filename
   * @returns MolecularStructure for rendering
   */
  public fromCP2K(content: string, filename?: string): MolecularStructure {
    const atoms = this.molecule3D.parseAtoms(content, 'CP2K');

    return {
      atoms,
      name: filename || 'CP2K Input',
    };
  }

  /**
   * Convert Quantum ESPRESSO input to MolecularStructure
   *
   * @param content - QE input file content
   * @param filename - Optional filename
   * @returns MolecularStructure for rendering
   */
  public fromQuantumEspresso(content: string, filename?: string): MolecularStructure {
    const atoms = this.molecule3D.parseAtoms(content, 'Quantum ESPRESSO');

    return {
      atoms,
      name: filename || 'QE Input',
    };
  }

  /**
   * Convert XYZ format to MolecularStructure
   *
   * @param content - XYZ file content
   * @param filename - Optional filename
   * @returns MolecularStructure for rendering
   */
  public fromXYZ(content: string, filename?: string): MolecularStructure {
    const lines = content.split('\n');
    const atoms: Atom[] = [];

    // XYZ format:
    // Line 1: Number of atoms
    // Line 2: Comment
    // Lines 3+: Element X Y Z

    if (lines.length < 3) {
      throw new Error('Invalid XYZ format: file too short');
    }

    const numAtoms = parseInt(lines[0].trim());
    if (isNaN(numAtoms)) {
      throw new Error('Invalid XYZ format: first line must be number of atoms');
    }

    for (let i = 2; i < Math.min(lines.length, numAtoms + 2); i++) {
      const line = lines[i].trim();
      if (!line || line.startsWith('#')) continue;

      const parts = line.split(/\s+/);
      if (parts.length >= 4) {
        atoms.push({
          element: parts[0],
          x: parseFloat(parts[1]),
          y: parseFloat(parts[2]),
          z: parseFloat(parts[3]),
        });
      }
    }

    return {
      atoms,
      name: filename || 'XYZ Structure',
    };
  }

  /**
   * Auto-detect format and convert to MolecularStructure
   *
   * @param content - File content
   * @param filename - Filename for format detection
   * @returns MolecularStructure for rendering
   */
  public autoConvert(content: string, filename: string): MolecularStructure {
    const ext = filename.split('.').pop()?.toLowerCase();

    // Use filename extension to determine format
    switch (ext) {
      case 'poscar':
      case 'contcar':
        return this.fromPOSCAR(content, filename);

      case 'gjf':
      case 'com':
        return this.fromGaussian(content, filename);

      case 'inp':
        // Could be ORCA or CP2K, try to detect from content
        if (content.includes('* xyz') || content.includes('* xyzfile')) {
          return this.fromORCA(content, filename);
        } else if (content.includes('&COORD')) {
          return this.fromCP2K(content, filename);
        }
        // Default to ORCA
        return this.fromORCA(content, filename);

      case 'xyz':
        return this.fromXYZ(content, filename);

      default:
        // Try to detect from content
        if (content.includes('&COORD')) {
          return this.fromCP2K(content, filename);
        } else if (content.includes('* xyz')) {
          return this.fromORCA(content, filename);
        } else if (content.includes('%chk')) {
          return this.fromGaussian(content, filename);
        } else {
          throw new Error(`Unknown file format: ${filename}`);
        }
    }
  }

  /**
   * Validate molecular structure
   *
   * @param structure - Structure to validate
   * @returns Validation result
   */
  public validateStructure(structure: MolecularStructure): {
    valid: boolean;
    errors: string[];
    warnings: string[];
  } {
    const errors: string[] = [];
    const warnings: string[] = [];

    // Check if atoms array exists
    if (!structure.atoms || structure.atoms.length === 0) {
      errors.push('No atoms found in structure');
    }

    // Check each atom
    structure.atoms?.forEach((atom, index) => {
      if (!atom.element) {
        errors.push(`Atom ${index}: missing element`);
      }

      if (
        isNaN(atom.x) ||
        isNaN(atom.y) ||
        isNaN(atom.z) ||
        !isFinite(atom.x) ||
        !isFinite(atom.y) ||
        !isFinite(atom.z)
      ) {
        errors.push(`Atom ${index}: invalid coordinates`);
      }
    });

    // Check for duplicate atoms (very close positions)
    const tolerance = 0.01; // Angstroms
    for (let i = 0; i < (structure.atoms?.length || 0); i++) {
      for (let j = i + 1; j < (structure.atoms?.length || 0); j++) {
        const atom1 = structure.atoms![i];
        const atom2 = structure.atoms![j];

        const dx = atom1.x - atom2.x;
        const dy = atom1.y - atom2.y;
        const dz = atom1.z - atom2.z;
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);

        if (dist < tolerance) {
          warnings.push(
            `Atoms ${i} and ${j} are very close (${dist.toFixed(3)} Å)`
          );
        }
      }
    }

    return {
      valid: errors.length === 0,
      errors,
      warnings,
    };
  }
}
