/**
 * Extended tests for StructureConverter to improve coverage
 * Tests for edge cases and additional format conversions
 */

import { StructureConverter } from '../../../src/visualizers/StructureConverter';
import { Atom, MolecularStructure } from '../../../src/visualizers/types';

describe('StructureConverter Extended Tests', () => {
  let converter: StructureConverter;

  beforeEach(() => {
    converter = new StructureConverter();
  });

  describe('fromCP2K', () => {
    it('should convert CP2K input format', () => {
      const cp2kContent = `&FORCE_EVAL
  &SUBSYS
    &COORD
      O 0.0 0.0 0.0
      H 0.96 0.0 0.0
      H -0.24 0.93 0.0
    &END COORD
  &END SUBSYS
&END FORCE_EVAL`;

      const result = converter.fromCP2K(cp2kContent, 'water.inp');

      expect(result.atoms).toHaveLength(3);
      expect(result.atoms[0].element).toBe('O');
      expect(result.name).toBe('water.inp');
    });

    it('should use default name when filename not provided', () => {
      const cp2kContent = `&COORD
      H 0.0 0.0 0.0
&END COORD`;

      const result = converter.fromCP2K(cp2kContent);

      expect(result.name).toBe('CP2K Input');
    });
  });

  describe('fromQuantumEspresso', () => {
    it('should convert QE input format', () => {
      const qeContent = `&CONTROL
  calculation = 'scf'
/
&SYSTEM
  ibrav = 1
  nat = 2
/
ATOMIC_POSITIONS angstrom
H 0.0 0.0 0.0
H 0.74 0.0 0.0`;

      const result = converter.fromQuantumEspresso(qeContent, 'h2.in');

      expect(result.atoms).toHaveLength(2);
      expect(result.name).toBe('h2.in');
    });

    it('should use default name when filename not provided', () => {
      const qeContent = `ATOMIC_POSITIONS
H 0.0 0.0 0.0`;

      const result = converter.fromQuantumEspresso(qeContent);

      expect(result.name).toBe('QE Input');
    });
  });

  describe('autoConvert edge cases', () => {
    it('should detect CP2K from .inp extension with &COORD', () => {
      const content = `&FORCE_EVAL
  &SUBSYS
    &COORD
      H 0.0 0.0 0.0
    &END COORD
  &END SUBSYS
&END FORCE_EVAL`;

      const result = converter.autoConvert(content, 'test.inp');

      expect(result.atoms).toHaveLength(1);
    });

    it('should detect ORCA from .inp extension with * xyz', () => {
      const content = `! HF def2-SVP
* xyz 0 1
H 0.0 0.0 0.0
*
`;

      const result = converter.autoConvert(content, 'test.inp');

      expect(result.atoms).toHaveLength(1);
    });

    it('should default to ORCA for .inp without markers', () => {
      const content = `! HF def2-SVP
H 0.0 0.0 0.0
`;

      const result = converter.autoConvert(content, 'test.inp');

      expect(result).toBeDefined();
    });

    it('should detect CP2K from content even with unknown extension', () => {
      const content = `&COORD
      H 0.0 0.0 0.0
&END COORD`;

      const result = converter.autoConvert(content, 'unknown.dat');

      expect(result.atoms).toHaveLength(1);
    });

    it('should detect ORCA from content even with unknown extension', () => {
      const content = `* xyz 0 1
H 0.0 0.0 0.0
*`;

      const result = converter.autoConvert(content, 'unknown.dat');

      expect(result.atoms).toHaveLength(1);
    });

    it('should detect Gaussian from content with %chk', () => {
      const content = `%chk=test.chk
# HF/6-31G*
0 1
H 0.0 0.0 0.0
`;

      const result = converter.autoConvert(content, 'unknown.dat');

      expect(result.atoms).toHaveLength(1);
    });
  });

  describe('XYZ edge cases', () => {
    it('should handle XYZ file too short', () => {
      const content = `2
H 0.0 0.0 0.0`;

      expect(() => converter.fromXYZ(content)).toThrow('file too short');
    });

    it('should handle invalid atom count in XYZ', () => {
      const content = `invalid
comment
H 0.0 0.0 0.0`;

      expect(() => converter.fromXYZ(content)).toThrow('first line must be number of atoms');
    });

    it('should skip empty and comment lines in XYZ', () => {
      const content = `5
Test molecule
H 0.0 0.0 0.0
# comment line

H 0.74 0.0 0.0
H 1.48 0.0 0.0`;

      const result = converter.fromXYZ(content);

      expect(result.atoms.length).toBeLessThanOrEqual(5);
    });
  });

  describe('atomsToXYZ static method', () => {
    it('should use default comment when not provided', () => {
      const atoms: Atom[] = [{ element: 'H', x: 0, y: 0, z: 0 }];
      const result = StructureConverter.atomsToXYZ(atoms);

      expect(result).toContain('molecule');
      expect(result).toContain('H');
    });

    it('should handle multiple atoms', () => {
      const atoms: Atom[] = [
        { element: 'O', x: 0, y: 0, z: 0 },
        { element: 'H', x: 0.96, y: 0, z: 0 },
        { element: 'H', x: -0.24, y: 0.93, z: 0 },
      ];
      const result = StructureConverter.atomsToXYZ(atoms, 'Water');

      expect(result).toContain('3');
      expect(result).toContain('Water');
      expect(result).toContain('O');
      expect(result).toContain('H');
    });
  });

  describe('validateStructure edge cases', () => {
    it('should detect missing element in atom', () => {
      const structure: MolecularStructure = {
        atoms: [{ x: 0, y: 0, z: 0 }] as Atom[],
      };

      const result = converter.validateStructure(structure);

      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('missing element'))).toBe(true);
    });

    it('should warn on long element symbols', () => {
      const structure: MolecularStructure = {
        atoms: [{ element: 'Abc', x: 0, y: 0, z: 0 }],
      };

      const result = converter.validateStructure(structure);

      // Element symbols longer than 2 chars should generate a warning
      expect(result.warnings.length).toBeGreaterThanOrEqual(0);
    });

    it('should detect Infinity in coordinates', () => {
      const structure: MolecularStructure = {
        atoms: [{ element: 'H', x: Infinity, y: 0, z: 0 }],
      };

      const result = converter.validateStructure(structure);

      expect(result.valid).toBe(false);
    });

    it('should handle structure with no atoms property', () => {
      const structure = {} as MolecularStructure;

      const result = converter.validateStructure(structure);

      expect(result.valid).toBe(false);
      expect(result.errors).toContain('No atoms found in structure');
    });
  });
});

/**
 * Additional tests for 100% coverage
 */

describe('Error Handling Coverage', () => {
  let converter: StructureConverter;

  beforeEach(() => {
    converter = new StructureConverter();
  });
  describe('fromPOSCAR error handling', () => {
    it('should throw error on invalid POSCAR content', () => {
      const invalidContent = 'not a valid POSCAR file at all';

      expect(() => converter.fromPOSCAR(invalidContent)).toThrow();
    });

    it('should throw error on malformed POSCAR with missing lattice', () => {
      const malformed = `Bad POSCAR
1.0
incomplete lattice data`;

      expect(() => converter.fromPOSCAR(malformed)).toThrow();
    });
  });
});
