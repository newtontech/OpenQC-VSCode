/**
 * Unit tests for OpenQCStructure validation.
 */

import { validateOpenQCStructure } from '../../../src/structures/validation';
import { OPENQC_STRUCTURE_SCHEMA_VERSION } from '../../../src/structures/OpenQCStructure';
import {
  WATER_STRUCTURE,
  SILICON_STRUCTURE,
} from '../../../src/structures/fixtures';

describe('validateOpenQCStructure', () => {
  // -----------------------------------------------------------------------
  // Valid structures
  // -----------------------------------------------------------------------

  describe('valid structures', () => {
    it('accepts a valid molecule structure (water)', () => {
      const result = validateOpenQCStructure(WATER_STRUCTURE);
      expect(result.valid).toBe(true);
      expect(result.errors).toEqual([]);
    });

    it('accepts a valid periodic structure (silicon)', () => {
      const result = validateOpenQCStructure(SILICON_STRUCTURE);
      expect(result.valid).toBe(true);
      expect(result.errors).toEqual([]);
    });

    it('accepts a minimal valid structure with just required fields', () => {
      const minimal = {
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      };
      const result = validateOpenQCStructure(minimal);
      expect(result.valid).toBe(true);
      expect(result.errors).toEqual([]);
    });
  });

  // -----------------------------------------------------------------------
  // Invalid schema version
  // -----------------------------------------------------------------------

  describe('schemaVersion validation', () => {
    it('rejects missing schemaVersion', () => {
      const result = validateOpenQCStructure({
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors).toContain('schemaVersion is required and must be a string');
    });

    it('rejects unknown schemaVersion', () => {
      const result = validateOpenQCStructure({
        schemaVersion: 'openqc.structure.v99',
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('Unsupported schemaVersion'))).toBe(true);
    });
  });

  // -----------------------------------------------------------------------
  // Invalid kind
  // -----------------------------------------------------------------------

  describe('kind validation', () => {
    it('rejects missing kind', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors).toContain('kind is required and must be a string');
    });

    it('rejects unknown kind', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'unknown_type',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('Invalid kind'))).toBe(true);
    });

    it('accepts all valid kinds', () => {
      const kinds = ['molecule', 'periodic', 'surface', 'trajectory', 'volumetric'];
      for (const kind of kinds) {
        const result = validateOpenQCStructure({
          schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
          kind,
          atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
        });
        expect(result.valid).toBe(true);
      }
    });
  });

  // -----------------------------------------------------------------------
  // Atom validation
  // -----------------------------------------------------------------------

  describe('atom validation', () => {
    it('rejects missing atoms array', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
      });
      expect(result.valid).toBe(false);
      expect(result.errors).toContain('atoms is required and must be an array');
    });

    it('rejects empty atoms array', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [],
      });
      expect(result.valid).toBe(false);
      expect(result.errors).toContain('atoms array must not be empty');
    });

    it('rejects atom with missing element', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ x: 0, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('missing or empty element'))).toBe(true);
    });

    it('rejects atom with empty element string', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: '  ', x: 0, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('missing or empty element'))).toBe(true);
    });

    it('rejects atom with NaN coordinates', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: NaN, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('x coordinate'))).toBe(true);
    });

    it('rejects atom with Infinity coordinates', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: Infinity, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('x coordinate'))).toBe(true);
    });

    it('rejects atom with invalid force tuple', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0, force: [1, 'bad'] as any }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('force must be a 3-tuple'))).toBe(true);
    });

    it('rejects atom with invalid selectiveDynamics', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0, selectiveDynamics: [true, false, 1] as any }],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('selectiveDynamics'))).toBe(true);
    });

    it('accepts atom with valid force tuple', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0, force: [0.1, 0.2, 0.3] }],
      });
      expect(result.valid).toBe(true);
    });

    it('accepts atom with valid selectiveDynamics', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0, selectiveDynamics: [true, false, true] }],
      });
      expect(result.valid).toBe(true);
    });
  });

  // -----------------------------------------------------------------------
  // Cell validation
  // -----------------------------------------------------------------------

  describe('cell validation', () => {
    it('rejects cell with missing lattice vectors', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'periodic',
        atoms: [{ element: 'Si', x: 0, y: 0, z: 0 }],
        cell: { pbc: [true, true, true] },
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('cell.a'))).toBe(true);
      expect(result.errors.some(e => e.includes('cell.b'))).toBe(true);
      expect(result.errors.some(e => e.includes('cell.c'))).toBe(true);
    });

    it('rejects cell with zero-length lattice vector', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'periodic',
        atoms: [{ element: 'Si', x: 0, y: 0, z: 0 }],
        cell: {
          a: [0, 0, 0],
          b: [1, 0, 0],
          c: [0, 1, 0],
          pbc: [true, true, true],
        },
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('cell.a has zero length'))).toBe(true);
    });

    it('rejects cell with invalid pbc', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'periodic',
        atoms: [{ element: 'Si', x: 0, y: 0, z: 0 }],
        cell: {
          a: [1, 0, 0],
          b: [0, 1, 0],
          c: [0, 0, 1],
          pbc: [1, 1, 1] as any,
        },
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('cell.pbc'))).toBe(true);
    });

    it('rejects cell with invalid coordinateMode', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'periodic',
        atoms: [{ element: 'Si', x: 0, y: 0, z: 0 }],
        cell: {
          a: [1, 0, 0],
          b: [0, 1, 0],
          c: [0, 0, 1],
          pbc: [true, true, true],
          coordinateMode: 'invalid',
        },
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('coordinateMode'))).toBe(true);
    });

    it('accepts molecule structure without cell', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
      });
      expect(result.valid).toBe(true);
    });
  });

  // -----------------------------------------------------------------------
  // Bond validation
  // -----------------------------------------------------------------------

  describe('bond validation', () => {
    it('rejects bond with non-number from/to', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }, { element: 'O', x: 1, y: 0, z: 0 }],
        bonds: [{ from: '0', to: 1 } as any],
      });
      expect(result.valid).toBe(false);
      expect(result.errors.some(e => e.includes('"from" and "to"'))).toBe(true);
    });
  });

  // -----------------------------------------------------------------------
  // Warnings
  // -----------------------------------------------------------------------

  describe('warnings', () => {
    it('warns about overlapping atoms by default', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'molecule',
        atoms: [
          { element: 'H', x: 0, y: 0, z: 0 },
          { element: 'H', x: 0.001, y: 0, z: 0 },
        ],
      });
      expect(result.valid).toBe(true);
      expect(result.warnings.some(w => w.includes('very close'))).toBe(true);
    });

    it('does not warn about overlapping atoms when checkOverlap is false', () => {
      const result = validateOpenQCStructure(
        {
          schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
          kind: 'molecule',
          atoms: [
            { element: 'H', x: 0, y: 0, z: 0 },
            { element: 'H', x: 0.001, y: 0, z: 0 },
          ],
        },
        { checkOverlap: false }
      );
      expect(result.valid).toBe(true);
      expect(result.warnings).toEqual([]);
    });

    it('warns about frame with no atoms', () => {
      const result = validateOpenQCStructure({
        schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
        kind: 'trajectory',
        atoms: [{ element: 'H', x: 0, y: 0, z: 0 }],
        frames: [{ atoms: [] }],
      });
      expect(result.valid).toBe(true);
      expect(result.warnings.some(w => w.includes('Frame 0'))).toBe(true);
    });
  });

  // -----------------------------------------------------------------------
  // Non-object input
  // -----------------------------------------------------------------------

  describe('non-object input', () => {
    it('rejects null', () => {
      const result = validateOpenQCStructure(null);
      expect(result.valid).toBe(false);
      expect(result.errors[0]).toContain('non-null object');
    });

    it('rejects a string', () => {
      const result = validateOpenQCStructure('not a structure');
      expect(result.valid).toBe(false);
    });

    it('rejects a number', () => {
      const result = validateOpenQCStructure(42);
      expect(result.valid).toBe(false);
    });
  });
});
