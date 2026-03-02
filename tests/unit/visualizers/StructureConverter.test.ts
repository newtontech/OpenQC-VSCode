/**
 * Tests for StructureConverter - Phase 2: 3D Visualization
 *
 * Tests the conversion of various quantum chemistry formats
 * to the MolecularStructure format used by ThreeJsRenderer.
 */

import { StructureConverter } from '../../../src/visualizers/StructureConverter';
import { Atom, MolecularStructure } from '../../../src/visualizers/types';

describe('StructureConverter', () => {
  let converter: StructureConverter;

  beforeEach(() => {
    converter = new StructureConverter();
  });

  describe('POSCAR conversion', () => {
    it('should convert a simple POSCAR file', () => {
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

      const result = converter.fromPOSCAR(poscarContent, 'test.poscar');

      expect(result).toBeDefined();
      expect(result.atoms).toHaveLength(3);
      expect(result.name).toBe('test.poscar');
      expect(result.atoms[0].element).toBe('H');
      expect(result.atoms[2].element).toBe('O');
    });

    it('should extract lattice vectors from POSCAR', () => {
      const poscarContent = `Simple cubic
1.0
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
H
1
Cartesian
0.0 0.0 0.0
`;

      const result = converter.fromPOSCAR(poscarContent);

      expect(result.lattice).toBeDefined();
      expect(result.lattice?.a).toEqual([5.0, 0.0, 0.0]);
      expect(result.lattice?.b).toEqual([0.0, 5.0, 0.0]);
      expect(result.lattice?.c).toEqual([0.0, 0.0, 5.0]);
    });

    it('should handle Direct coordinate mode', () => {
      const poscarContent = `Test
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
H
1
Direct
0.5 0.5 0.5
`;

      const result = converter.fromPOSCAR(poscarContent);

      expect(result.atoms).toHaveLength(1);
      expect(result.atoms[0].element).toBe('H');
    });

    it('uses default comment when not provided', () => {
      const atoms = [{ element: 'C', x: 0, y: 0, z: 0 }];
      // Call without comment parameter to test default value branch
      const result = StructureConverter.atomsToXYZ(atoms);

      expect(result).toContain('molecule');
      expect(result).toContain('C');
    });
  });

  describe('XYZ conversion', () => {
    it('should convert a simple XYZ file', () => {
      const xyzContent = `3
Water molecule
O 0.0 0.0 0.0
H 0.96 0.0 0.0
H -0.24 0.93 0.0
`;

      const result = converter.fromXYZ(xyzContent, 'water.xyz');

      expect(result.atoms).toHaveLength(3);
      expect(result.atoms[0].element).toBe('O');
      expect(result.atoms[1].element).toBe('H');
      expect(result.atoms[2].element).toBe('H');
      expect(result.name).toBe('water.xyz');
    });

    it('should handle XYZ with comment lines', () => {
      const xyzContent = `2
# H2 molecule
H 0.0 0.0 0.0
H 0.0 0.0 0.74
`;

      const result = converter.fromXYZ(xyzContent);

      expect(result.atoms).toHaveLength(2);
    });

    it('should throw error for invalid XYZ format', () => {
      const invalidContent = `invalid
content
here`;

      expect(() => converter.fromXYZ(invalidContent)).toThrow();
    });
  });

  describe('Gaussian conversion', () => {
    it('should convert Gaussian input format', () => {
      const gaussianContent = `%chk=test.chk
# HF/6-31G*

Water

0 1
O 0.0 0.0 0.0
H 0.96 0.0 0.0
H -0.24 0.93 0.0
`;

      const result = converter.fromGaussian(gaussianContent, 'water.gjf');

      expect(result.atoms).toHaveLength(3);
      expect(result.atoms[0].element).toBe('O');
    });
  });

  describe('ORCA conversion', () => {
    it('should convert ORCA input format', () => {
      const orcaContent = `! HF def2-SVP
* xyz 0 1
O 0.0 0.0 0.0
H 0.96 0.0 0.0
H -0.24 0.93 0.0
*
`;

      const result = converter.fromORCA(orcaContent, 'water.inp');

      expect(result.atoms).toHaveLength(3);
      expect(result.atoms[0].element).toBe('O');
    });
  });

  describe('Structure validation', () => {
    it('should validate a correct structure', () => {
      const structure: MolecularStructure = {
        atoms: [
          { element: 'H', x: 0, y: 0, z: 0 },
          { element: 'O', x: 1, y: 0, z: 0 },
        ],
      };

      const result = converter.validateStructure(structure);

      expect(result.valid).toBe(true);
      expect(result.errors).toHaveLength(0);
    });

    it('should detect empty atoms array', () => {
      const structure: MolecularStructure = {
        atoms: [],
      };

      const result = converter.validateStructure(structure);

      expect(result.valid).toBe(false);
      expect(result.errors).toContain('No atoms found in structure');
    });

    it('should detect invalid coordinates', () => {
      const structure: MolecularStructure = {
        atoms: [{ element: 'H', x: NaN, y: 0, z: 0 }],
      };

      const result = converter.validateStructure(structure);

      expect(result.valid).toBe(false);
      expect(result.errors.length).toBeGreaterThan(0);
    });

    it('should warn about duplicate atoms', () => {
      const structure: MolecularStructure = {
        atoms: [
          { element: 'H', x: 0, y: 0, z: 0 },
          { element: 'H', x: 0.001, y: 0, z: 0 },
        ],
      };

      const result = converter.validateStructure(structure);

      expect(result.warnings.length).toBeGreaterThan(0);
      expect(result.warnings[0]).toContain('very close');
    });
  });

  describe('Auto-conversion', () => {
    it('should detect POSCAR format from filename', () => {
      const poscarContent = `Test
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
H
1
Cartesian
0.0 0.0 0.0
`;

      const result = converter.autoConvert(poscarContent, 'POSCAR');

      expect(result.atoms).toHaveLength(1);
    });

    it('should detect XYZ format from filename', () => {
      const xyzContent = `1
Atom
H 0.0 0.0 0.0
`;

      const result = converter.autoConvert(xyzContent, 'test.xyz');

      expect(result.atoms).toHaveLength(1);
    });

    it('should detect Gaussian format from filename', () => {
      const gaussianContent = `# HF/6-31G*

0 1
H 0.0 0.0 0.0
`;

      const result = converter.autoConvert(gaussianContent, 'test.gjf');

      expect(result.atoms).toHaveLength(1);
    });

    it('should throw error for unknown format', () => {
      expect(() => converter.autoConvert('random content', 'unknown.dat')).toThrow();
    });
  });
});
