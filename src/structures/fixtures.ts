/**
 * Reusable test fixtures for OpenQCStructure.
 *
 * Provides sample structures that tests and development tools can import
 * instead of duplicating inline data.
 *
 * @module structures/fixtures
 */

import type { OpenQCStructure, OpenQCAtom, OpenQCCell } from './OpenQCStructure';
import { OPENQC_STRUCTURE_SCHEMA_VERSION } from './OpenQCStructure';

// ---------------------------------------------------------------------------
// Water molecule (molecule-only, no cell)
// ---------------------------------------------------------------------------

export const WATER_ATOMS: OpenQCAtom[] = [
  { element: 'O', x: 0.0, y: 0.0, z: 0.1173 },
  { element: 'H', x: 0.0, y: 0.7572, z: -0.4692 },
  { element: 'H', x: 0.0, y: -0.7572, z: -0.4692 },
];

export const WATER_STRUCTURE: OpenQCStructure = {
  schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
  kind: 'molecule',
  name: 'Water',
  atoms: WATER_ATOMS,
  metadata: {
    charge: 0,
    multiplicity: 1,
    source: { software: 'test', parser: 'fixture' },
    provenance: {
      createdAt: '2025-01-01T00:00:00.000Z',
      warnings: [],
    },
  },
};

// ---------------------------------------------------------------------------
// Silicon crystal (periodic, with cell)
// ---------------------------------------------------------------------------

export const SILICON_CELL: OpenQCCell = {
  a: [5.43, 0.0, 0.0],
  b: [0.0, 5.43, 0.0],
  c: [0.0, 0.0, 5.43],
  pbc: [true, true, true],
  coordinateMode: 'fractional',
};

export const SILICON_ATOMS: OpenQCAtom[] = [
  { element: 'Si', x: 0.0, y: 0.0, z: 0.0 },
  { element: 'Si', x: 0.25, y: 0.25, z: 0.25 },
];

export const SILICON_STRUCTURE: OpenQCStructure = {
  schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
  kind: 'periodic',
  name: 'Silicon',
  atoms: SILICON_ATOMS,
  cell: SILICON_CELL,
  metadata: {
    source: { software: 'vasp', parser: 'poscar' },
    provenance: {
      createdAt: '2025-01-01T00:00:00.000Z',
      warnings: [],
    },
  },
};

// ---------------------------------------------------------------------------
// POSCAR-like content string (for testing poscarToOpenQCStructure)
// ---------------------------------------------------------------------------

export const SI_POSCAR_CONTENT = `Silicon
5.43
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
`;

// ---------------------------------------------------------------------------
// XYZ content strings
// ---------------------------------------------------------------------------

export const WATER_XYZ_CONTENT = `3
Water molecule
O   0.000000  0.000000  0.117300
H   0.000000  0.757200 -0.469200
H   0.000000 -0.757200 -0.469200
`;

export const METHANE_XYZ_CONTENT = `5
Methane
C   0.000000  0.000000  0.000000
H   0.627600  0.627600  0.627600
H  -0.627600 -0.627600  0.627600
H  -0.627600  0.627600 -0.627600
H   0.627600 -0.627600 -0.627600
`;
