/**
 * Scientific regression tests comparing parsed outputs against golden fixtures.
 *
 * These tests ensure that parser changes do not silently alter coordinates,
 * atom counts, element ordering, or lattice vectors.
 *
 * Issue #81: Golden fixtures, CI checks, and dependency/license review.
 *
 * @module tests/unit/fixtures/scientificRegression.test
 */

import * as fs from 'fs';
import * as path from 'path';
import { xyzToOpenQCStructure, poscarToOpenQCStructure } from '../../../src/structures/converters';
import type { OpenQCStructure } from '../../../src/structures/OpenQCStructure';
import { validateOpenQCStructure } from '../../../src/structures/validation';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const FIXTURES_ROOT = path.resolve(__dirname, '../../fixtures/structures');
const GOLDEN_ROOT = path.resolve(__dirname, '../../golden/structures');

interface GoldenStructure {
  schemaVersion: string;
  kind: string;
  name?: string;
  atomCount: number;
  elements: string[];
  cell?: {
    a: [number, number, number];
    b: [number, number, number];
    c: [number, number, number];
  };
  atoms?: Array<{ element: string; x: number; y: number; z: number }>;
}

function readFixture(relPath: string): string {
  const fullPath = path.join(FIXTURES_ROOT, relPath);
  return fs.readFileSync(fullPath, 'utf-8');
}

function readGolden(name: string): GoldenStructure {
  const fullPath = path.join(GOLDEN_ROOT, name);
  return JSON.parse(fs.readFileSync(fullPath, 'utf-8'));
}

/** Compare two numbers within tolerance (Angstroms). */
function approxEqual(a: number, b: number, tolerance: number = 1e-4): boolean {
  return Math.abs(a - b) < tolerance;
}

// ---------------------------------------------------------------------------
// XYZ Tests
// ---------------------------------------------------------------------------

describe('XYZ golden regression tests', () => {
  it('water.xyz matches golden', () => {
    const content = readFixture('xyz/water.xyz');
    const structure = xyzToOpenQCStructure(content, { name: 'Water molecule' });
    const golden = readGolden('water-xyz.json');

    expect(structure.atoms.length).toBe(golden.atomCount);
    expect(structure.kind).toBe(golden.kind);

    for (let i = 0; i < structure.atoms.length; i++) {
      expect(structure.atoms[i].element).toBe(golden.elements[i]);
      if (golden.atoms) {
        expect(approxEqual(structure.atoms[i].x, golden.atoms[i].x)).toBe(true);
        expect(approxEqual(structure.atoms[i].y, golden.atoms[i].y)).toBe(true);
        expect(approxEqual(structure.atoms[i].z, golden.atoms[i].z)).toBe(true);
      }
    }

    // Validate DTO
    const validation = validateOpenQCStructure(structure);
    expect(validation.valid).toBe(true);
  });

  it('methane.xyz matches golden', () => {
    const content = readFixture('xyz/methane.xyz');
    const structure = xyzToOpenQCStructure(content, { name: 'Methane' });
    const golden = readGolden('methane-xyz.json');

    expect(structure.atoms.length).toBe(golden.atomCount);
    expect(structure.atoms[0].element).toBe('C');
    expect(structure.kind).toBe('molecule');

    const validation = validateOpenQCStructure(structure);
    expect(validation.valid).toBe(true);
  });

  it('benzene.xyz matches golden atom count and elements', () => {
    const content = readFixture('xyz/benzene.xyz');
    const structure = xyzToOpenQCStructure(content);
    const golden = readGolden('benzene-xyz.json');

    expect(structure.atoms.length).toBe(golden.atomCount);
    expect(structure.atoms.map(a => a.element)).toEqual(golden.elements);
  });
});

// ---------------------------------------------------------------------------
// POSCAR Tests
// ---------------------------------------------------------------------------

describe('POSCAR golden regression tests', () => {
  it('Si_direct.vasp matches golden', () => {
    const content = readFixture('poscar/Si_direct.vasp');
    const structure = poscarToOpenQCStructure(content, 'Si_direct.vasp');
    const golden = readGolden('si-poscar-direct.json');

    expect(structure.atoms.length).toBe(golden.atomCount);
    expect(structure.kind).toBe('periodic');
    expect(structure.cell).toBeDefined();
    expect(structure.cell!.coordinateMode).toBe('fractional');

    // Check lattice vectors
    if (golden.cell) {
      for (const key of ['a', 'b', 'c'] as const) {
        for (let i = 0; i < 3; i++) {
          expect(approxEqual(structure.cell![key][i], golden.cell[key][i], 0.01)).toBe(true);
        }
      }
    }

    // Check atoms
    for (let i = 0; i < structure.atoms.length; i++) {
      expect(structure.atoms[i].element).toBe(golden.elements[i]);
    }

    const validation = validateOpenQCStructure(structure);
    expect(validation.valid).toBe(true);
  });

  it('NaCl_cart.vasp matches golden', () => {
    const content = readFixture('poscar/NaCl_cart.vasp');
    const structure = poscarToOpenQCStructure(content, 'NaCl_cart.vasp');
    const golden = readGolden('nacl-poscar-cart.json');

    expect(structure.atoms.length).toBe(golden.atomCount);
    expect(structure.kind).toBe('periodic');
    expect(structure.cell).toBeDefined();
    expect(structure.cell!.coordinateMode).toBe('cartesian');

    // Check atoms
    for (let i = 0; i < structure.atoms.length; i++) {
      expect(structure.atoms[i].element).toBe(golden.elements[i]);
      if (golden.atoms) {
        expect(approxEqual(structure.atoms[i].x, golden.atoms[i].x, 0.01)).toBe(true);
        expect(approxEqual(structure.atoms[i].y, golden.atoms[i].y, 0.01)).toBe(true);
        expect(approxEqual(structure.atoms[i].z, golden.atoms[i].z, 0.01)).toBe(true);
      }
    }
  });

  it('Fe_selective.vasp parses selective dynamics flags', () => {
    const content = readFixture('poscar/Fe_selective.vasp');
    const structure = poscarToOpenQCStructure(content, 'Fe_selective.vasp');

    expect(structure.atoms.length).toBe(2);
    expect(structure.atoms[0].selectiveDynamics).toEqual([true, true, true]);
    expect(structure.atoms[1].selectiveDynamics).toEqual([false, false, true]);
  });
});

// ---------------------------------------------------------------------------
// DTO Validation Tests
// ---------------------------------------------------------------------------

describe('Structure DTO validation with fixtures', () => {
  it('water XYZ produces valid DTO with correct metadata', () => {
    const content = readFixture('xyz/water.xyz');
    const structure = xyzToOpenQCStructure(content);

    expect(structure.schemaVersion).toBe('openqc.structure.v1');
    expect(structure.metadata?.source?.parser).toBe('xyz');
    expect(structure.metadata?.provenance?.createdAt).toBeDefined();
    expect(structure.metadata?.provenance?.warnings).toEqual([]);
  });

  it('POSCAR produces valid periodic DTO', () => {
    const content = readFixture('poscar/Si_direct.vasp');
    const structure = poscarToOpenQCStructure(content);

    expect(structure.schemaVersion).toBe('openqc.structure.v1');
    expect(structure.kind).toBe('periodic');
    expect(structure.metadata?.source?.software).toBe('vasp');
    expect(structure.metadata?.source?.parser).toBe('poscar');
    expect(structure.cell?.pbc).toEqual([true, true, true]);
  });
});

// ---------------------------------------------------------------------------
// Fixture Completeness Checks
// ---------------------------------------------------------------------------

describe('Fixture directory completeness', () => {
  const expectedFixtureDirs = ['poscar', 'cif', 'xyz', 'qe', 'cp2k', 'gaussian', 'orca'];

  for (const dir of expectedFixtureDirs) {
    it(`${dir}/ directory exists and has at least one fixture`, () => {
      const dirPath = path.join(FIXTURES_ROOT, dir);
      expect(fs.existsSync(dirPath)).toBe(true);

      const files = fs.readdirSync(dirPath).filter(f => !f.startsWith('.'));
      expect(files.length).toBeGreaterThanOrEqual(1);
    });
  }

  it('golden directory exists with at least 4 golden files', () => {
    expect(fs.existsSync(GOLDEN_ROOT)).toBe(true);
    const files = fs.readdirSync(GOLDEN_ROOT).filter(f => f.endsWith('.json'));
    expect(files.length).toBeGreaterThanOrEqual(4);
  });
});
