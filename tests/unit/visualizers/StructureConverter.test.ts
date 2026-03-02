import { StructureConverter } from '../../../src/visualizers/StructureConverter';

describe('StructureConverter', () => {
  describe('atomsToXYZ', () => {
    it('converts atoms array to XYZ format string', () => {
      const atoms = [
        { elem: 'H', x: 0.0, y: 0.0, z: 0.0 },
        { elem: 'O', x: 0.0, y: 0.0, z: 1.0 },
      ];

      const result = StructureConverter.atomsToXYZ(atoms, 'water');

      expect(result).toContain('2');
      expect(result).toContain('water');
      expect(result).toMatch(/H\s+0\.0+\s+0\.0+\s+0\.0+/);
      expect(result).toMatch(/O\s+0\.0+\s+0\.0+\s+1\.0+/);
    });

    it('handles empty atoms array', () => {
      const result = StructureConverter.atomsToXYZ([], 'empty');

      expect(result).toContain('0');
      expect(result).toContain('empty');
    });

    it('formats coordinates with consistent precision', () => {
      const atoms = [{ elem: 'C', x: 0.123456789, y: 1.23456789, z: 2.3456789 }];

      const result = StructureConverter.atomsToXYZ(atoms, 'carbon');

      // Should have consistent decimal places (6 is typical for XYZ)
      const lines = result.split('\n');
      const coordLine = lines[2]; // First atom line
      const parts = coordLine.trim().split(/\s+/);

      expect(parseFloat(parts[1])).toBeCloseTo(0.123457, 5);
      expect(parseFloat(parts[2])).toBeCloseTo(1.234568, 5);
      expect(parseFloat(parts[3])).toBeCloseTo(2.345679, 5);
    });

    it('handles special characters in comment', () => {
      const atoms = [{ elem: 'Fe', x: 0, y: 0, z: 0 }];
      const result = StructureConverter.atomsToXYZ(atoms, 'Iron oxide (Fe2O3)');

      expect(result).toContain('Iron oxide (Fe2O3)');
    });
  });

  describe('atomsToJSON', () => {
    it('converts atoms array to NGL-compatible JSON', () => {
      const atoms = [
        { elem: 'H', x: 0, y: 0, z: 0 },
        { elem: 'O', x: 0, y: 0, z: 1 },
      ];

      const result = StructureConverter.atomsToJSON(atoms);
      const parsed = JSON.parse(result);

      expect(parsed).toHaveProperty('atoms');
      expect(parsed.atoms).toHaveLength(2);
      expect(parsed.atoms[0]).toEqual({ elem: 'H', x: 0, y: 0, z: 0 });
      expect(parsed.atoms[1]).toEqual({ elem: 'O', x: 0, y: 0, z: 1 });
    });

    it('includes unit cell data when provided', () => {
      const atoms = [{ elem: 'Si', x: 0, y: 0, z: 0 }];
      const unitCell = {
        a: 5.43,
        b: 5.43,
        c: 5.43,
        alpha: 90,
        beta: 90,
        gamma: 90,
      };

      const result = StructureConverter.atomsToJSON(atoms, unitCell);
      const parsed = JSON.parse(result);

      expect(parsed).toHaveProperty('unitCell');
      expect(parsed.unitCell).toEqual(unitCell);
    });
  });
});
