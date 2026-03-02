/**
 * Integration tests for Three.js visualization
 *
 * Tests the complete pipeline from file parsing to visualization
 */

import { StructureConverter } from '../../../src/visualizers/StructureConverter';
import { VASPParser } from '../../../src/parsers/VASPParser';
import fs from 'fs';
import path from 'path';

describe('Three.js Visualization Integration', () => {
  let converter: StructureConverter;

  beforeEach(() => {
    converter = new StructureConverter();
  });

  describe('VASP POSCAR to Visualization Pipeline', () => {
    it('should convert POSCAR file to renderable structure', () => {
      const poscarContent = `H2O molecule
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
H O
2 1
Cartesian
0.0 0.0 0.0
0.0 0.0 0.74
0.5 0.5 0.5
`;

      const structure = converter.fromPOSCAR(poscarContent, 'H2O.poscar');

      // Verify structure was converted
      expect(structure).toBeDefined();
      expect(structure.atoms).toHaveLength(3);
      expect(structure.name).toBe('H2O.poscar');

      // Verify atom data is correct
      const oxygenAtom = structure.atoms.find((a) => a.element === 'O');
      expect(oxygenAtom).toBeDefined();
      expect(oxygenAtom?.x).toBe(0.5);
      expect(oxygenAtom?.y).toBe(0.5);
      expect(oxygenAtom?.z).toBe(0.5);

      // Verify lattice vectors
      expect(structure.lattice).toBeDefined();
      expect(structure.lattice?.a).toEqual([10.0, 0.0, 0.0]);
      expect(structure.lattice?.b).toEqual([0.0, 10.0, 0.0]);
      expect(structure.lattice?.c).toEqual([0.0, 0.0, 10.0]);

      // Validate the structure
      const validation = converter.validateStructure(structure);
      expect(validation.valid).toBe(true);
      expect(validation.errors).toHaveLength(0);
    });

    it('should handle complex multi-element structures', () => {
      const poscarContent = `SiO2 quartz structure
1.0
4.913 0.0 0.0
-2.4565 4.2548 0.0
0.0 0.0 5.4052
Si O
3 6
Cartesian
0.0 0.0 0.0
2.4565 1.4183 1.8017
1.2282 2.8365 3.6035
1.2282 0.7091 0.9009
3.6848 2.1274 2.7026
-1.2282 2.1274 0.9009
1.2282 -0.7091 2.7026
3.6848 0.7091 3.6035
-1.2282 0.7091 1.8017
`;

      const structure = converter.fromPOSCAR(poscarContent);

      expect(structure.atoms).toHaveLength(9);
      expect(structure.atoms.filter((a) => a.element === 'Si')).toHaveLength(3);
      expect(structure.atoms.filter((a) => a.element === 'O')).toHaveLength(6);

      const validation = converter.validateStructure(structure);
      expect(validation.valid).toBe(true);
    });
  });

  describe('Format Detection and Auto-Conversion', () => {
    it('should auto-detect and convert POSCAR format', () => {
      const poscarContent = `Simple structure
1.0
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
C
1
Cartesian
0.0 0.0 0.0
`;

      const structure = converter.autoConvert(poscarContent, 'POSCAR');

      expect(structure.atoms).toHaveLength(1);
      expect(structure.atoms[0].element).toBe('C');
    });

    it('should auto-detect and convert XYZ format', () => {
      const xyzContent = `4
Methane
C 0.0 0.0 0.0
H 0.63 0.63 0.63
H -0.63 -0.63 0.63
H -0.63 0.63 -0.63
H 0.63 -0.63 -0.63
`;

      const structure = converter.autoConvert(xyzContent, 'methane.xyz');

      expect(structure.atoms).toHaveLength(5);
      expect(structure.atoms[0].element).toBe('C');
      expect(structure.atoms.filter((a) => a.element === 'H')).toHaveLength(4);
    });

    it('should auto-detect and convert Gaussian format', () => {
      const gaussianContent = `%chk=methane.chk
# B3LYP/6-31G*

Methane

0 1
C 0.0 0.0 0.0
H 0.63 0.63 0.63
H -0.63 -0.63 0.63
H -0.63 0.63 -0.63
H 0.63 -0.63 -0.63
`;

      const structure = converter.autoConvert(gaussianContent, 'methane.gjf');

      expect(structure.atoms).toHaveLength(5);
    });
  });

  describe('Structure Validation', () => {
    it('should detect and report structure issues', () => {
      // Create a structure with problematic atoms
      const poscarContent = `Problematic structure
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
2
Cartesian
0.0 0.0 0.0
0.001 0.001 0.001
`;

      const structure = converter.fromPOSCAR(poscarContent);
      const validation = converter.validateStructure(structure);

      // Structure should be valid (no critical errors)
      expect(validation.valid).toBe(true);

      // But should have warnings about duplicate atoms
      expect(validation.warnings.length).toBeGreaterThan(0);
      expect(validation.warnings[0]).toContain('very close');
    });

    it('should reject structures with invalid coordinates', () => {
      const poscarContent = `Invalid structure
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
1
Cartesian
invalid 0.0 0.0
`;

      const structure = converter.fromPOSCAR(poscarContent);
      const validation = converter.validateStructure(structure);

      // Should have errors about invalid coordinates
      expect(validation.errors.length).toBeGreaterThan(0);
    });
  });

  describe('VASP Parser Integration', () => {
    it('should work with VASPParser for coordinate extraction', () => {
      const poscarContent = `Test structure
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
H O
2 1
Cartesian
0.0 0.0 0.0
0.0 0.0 0.74
0.5 0.5 0.5
`;

      const parser = new VASPParser(poscarContent, 'POSCAR');

      // Test coordinate extraction
      const coordinates = parser.getCoordinates();
      expect(coordinates).toHaveLength(3);
      expect(coordinates[0]).toEqual([0, 0, 0]);
      expect(coordinates[1]).toEqual([0, 0, 0.74]);

      // Test lattice vector extraction
      const latticeVectors = parser.getLatticeVectors();
      expect(latticeVectors).toHaveLength(3);
      expect(latticeVectors[0]).toEqual([10.0, 0.0, 0.0]);

      // Test atom type extraction
      const atomTypes = parser.getAtomTypes();
      expect(atomTypes).toEqual(['H', 'O']);

      // Test atom count extraction
      const atomCounts = parser.getAtomCounts();
      expect(atomCounts).toEqual([2, 1]);
    });
  });
});
