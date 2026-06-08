/**
 * Conversion utilities to/from OpenQCStructure.
 *
 * Provides:
 * - `openQCStructureToXYZ`  -- backward-compatible XYZ string output.
 * - `xyzToOpenQCStructure`  -- parse XYZ text into the canonical DTO.
 * - `molecularStructureToOpenQCStructure` -- adapt existing MolecularStructure.
 *
 * @module structures/converters
 */

import type {
  OpenQCStructure,
  OpenQCAtom,
  OpenQCCell,
  OpenQCStructureKind,
} from './OpenQCStructure';
import { OPENQC_STRUCTURE_SCHEMA_VERSION } from './OpenQCStructure';
import type { Atom, MolecularStructure, LatticeVectors } from '../visualizers/types';
import { validateOpenQCStructure } from './validation';

// ---------------------------------------------------------------------------
// XYZ <-> OpenQCStructure
// ---------------------------------------------------------------------------

/**
 * Convert an OpenQCStructure to an XYZ-format string.
 *
 * This is the backward-compatible path used by the existing 3Dmol/NGL
 * viewers that expect an XYZ string.
 *
 * @param structure - The canonical structure DTO.
 * @returns XYZ string (atom-count header, comment line, then atom lines).
 */
export function openQCStructureToXYZ(structure: OpenQCStructure): string {
  const atoms = structure.atoms;
  const name = structure.name ?? 'OpenQCStructure';
  const lines: string[] = [String(atoms.length), name];

  for (const atom of atoms) {
    const element = atom.element.padEnd(2);
    const x = atom.x.toFixed(6);
    const y = atom.y.toFixed(6);
    const z = atom.z.toFixed(6);
    lines.push(`${element} ${x} ${y} ${z}`);
  }

  return lines.join('\n') + '\n';
}

/**
 * Parse an XYZ-format string into an OpenQCStructure.
 *
 * @param content  - XYZ file content.
 * @param options  - Optional overrides for kind, name, and source metadata.
 * @returns An OpenQCStructure of kind 'molecule'.
 * @throws Error if the XYZ content is malformed.
 */
export function xyzToOpenQCStructure(
  content: string,
  options?: {
    kind?: OpenQCStructureKind;
    name?: string;
    sourceSoftware?: string;
  }
): OpenQCStructure {
  const lines = content.split('\n');

  if (lines.length < 3) {
    throw new Error('Invalid XYZ format: file too short');
  }

  const numAtoms = parseInt(lines[0].trim(), 10);
  if (isNaN(numAtoms) || numAtoms <= 0) {
    throw new Error('Invalid XYZ format: first line must be a positive integer (atom count)');
  }

  const comment = lines[1].trim();
  const atoms: OpenQCAtom[] = [];

  for (let i = 2; i < Math.min(lines.length, numAtoms + 2); i++) {
    const line = lines[i].trim();
    if (!line || line.startsWith('#')) {
      continue;
    }

    const parts = line.split(/\s+/);
    if (parts.length < 4) {
      continue;
    }

    const element = parts[0];
    const x = parseFloat(parts[1]);
    const y = parseFloat(parts[2]);
    const z = parseFloat(parts[3]);

    if (isNaN(x) || isNaN(y) || isNaN(z)) {
      continue;
    }

    atoms.push({ element, x, y, z });
  }

  if (atoms.length === 0) {
    throw new Error('Invalid XYZ format: no valid atom lines found');
  }

  const structure: OpenQCStructure = {
    schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
    kind: options?.kind ?? 'molecule',
    name: options?.name ?? (comment || 'XYZ Structure'),
    atoms,
    metadata: {
      source: {
        software: options?.sourceSoftware,
        parser: 'xyz',
      },
      provenance: {
        createdAt: new Date().toISOString(),
        warnings: [],
      },
    },
  };

  return structure;
}

// ---------------------------------------------------------------------------
// MolecularStructure (existing visualizer type) <-> OpenQCStructure
// ---------------------------------------------------------------------------

function latticeVectorsToCell(lattice: LatticeVectors): OpenQCCell {
  return {
    a: lattice.a,
    b: lattice.b,
    c: lattice.c,
    pbc: [true, true, true],
  };
}

function cellToLatticeVectors(cell: OpenQCCell): LatticeVectors {
  return {
    a: cell.a,
    b: cell.b,
    c: cell.c,
  };
}

/**
 * Convert the existing `MolecularStructure` (from `visualizers/types`) into
 * the canonical `OpenQCStructure`.
 *
 * @param ms      - The legacy molecular structure.
 * @param options - Optional overrides for kind and source metadata.
 * @returns An OpenQCStructure.
 */
export function molecularStructureToOpenQCStructure(
  ms: MolecularStructure,
  options?: {
    kind?: OpenQCStructureKind;
    sourceSoftware?: string;
    sourceParser?: string;
  }
): OpenQCStructure {
  const hasLattice = ms.lattice !== undefined;
  const kind: OpenQCStructureKind = options?.kind ?? (hasLattice ? 'periodic' : 'molecule');

  const atoms: OpenQCAtom[] = ms.atoms.map(a => ({
    element: a.element,
    x: a.x,
    y: a.y,
    z: a.z,
  }));

  const structure: OpenQCStructure = {
    schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
    kind,
    name: ms.name,
    atoms,
    metadata: {
      source: {
        software: options?.sourceSoftware,
        parser: options?.sourceParser,
      },
      provenance: {
        createdAt: new Date().toISOString(),
        warnings: [],
      },
    },
  };

  if (hasLattice && ms.lattice) {
    structure.cell = latticeVectorsToCell(ms.lattice);
  }

  if (ms.bonds && ms.bonds.length > 0) {
    structure.bonds = ms.bonds.map(b => ({
      from: b.atomIndex1,
      to: b.atomIndex2,
      order: b.order,
    }));
  }

  return structure;
}

/**
 * Convert an OpenQCStructure back to the legacy `MolecularStructure`.
 *
 * This is the compatibility adapter that existing viewer code can use
 * while it has not yet been migrated to the new DTO.
 *
 * @param structure - The canonical structure DTO.
 * @returns A MolecularStructure compatible with the existing visualizers.
 */
export function openQCStructureToMolecularStructure(
  structure: OpenQCStructure
): MolecularStructure {
  const atoms: Atom[] = structure.atoms.map(a => ({
    element: a.element,
    x: a.x,
    y: a.y,
    z: a.z,
  }));

  const ms: MolecularStructure = {
    atoms,
    name: structure.name,
  };

  if (structure.cell) {
    ms.lattice = cellToLatticeVectors(structure.cell);
  }

  if (structure.bonds && structure.bonds.length > 0) {
    ms.bonds = structure.bonds.map(b => ({
      atomIndex1: b.from,
      atomIndex2: b.to,
      length: 0,
      order: b.order,
    }));
  }

  return ms;
}

// ---------------------------------------------------------------------------
// POSCAR <-> OpenQCStructure
// ---------------------------------------------------------------------------

/**
 * Convert a VASP POSCAR content string into an OpenQCStructure.
 *
 * @param content - POSCAR file content.
 * @param filename - Optional filename for metadata.
 * @returns An OpenQCStructure of kind 'periodic' (or 'molecule' if no lattice).
 * @throws Error if the POSCAR content cannot be parsed.
 */
export function poscarToOpenQCStructure(content: string, filename?: string): OpenQCStructure {
  const lines = content.split('\n');

  if (lines.length < 8) {
    throw new Error('Invalid POSCAR format: too few lines');
  }

  const comment = lines[0].trim();
  const scale = parseFloat(lines[1].trim());
  const scaleF = isNaN(scale) ? 1.0 : scale;

  // Lattice vectors
  const parseVec = (line: string): [number, number, number] => {
    const parts = line.trim().split(/\s+/);
    return [
      parseFloat(parts[0]) * scaleF,
      parseFloat(parts[1]) * scaleF,
      parseFloat(parts[2]) * scaleF,
    ];
  };

  const a = parseVec(lines[2]);
  const b = parseVec(lines[3]);
  const c = parseVec(lines[4]);

  // Atom types and counts
  const typeLine = lines[5].trim().split(/\s+/);
  const countLine = lines[6].trim().split(/\s+/).map(Number);

  // Determine if types are on line 5 or not (older POSCAR format)
  const hasElementNames = typeLine.every(t => /^[A-Za-z]+$/.test(t));
  const elementNames = hasElementNames ? typeLine : countLine.map((_, i) => `El${i + 1}`);
  const atomCounts = hasElementNames ? countLine : typeLine.map(Number);

  // Coordinate mode
  const modeLine = lines[7].trim().toLowerCase();
  const isSelective = modeLine.startsWith('s');
  const isDirect =
    modeLine.startsWith('d') || (!isSelective && !modeLine.startsWith('c') && modeLine !== '');
  const coordLine = isSelective || isDirect || modeLine.startsWith('c') ? 8 : 7;
  const coordinateMode = isDirect ? 'fractional' : 'cartesian';

  const cell: OpenQCCell = {
    a,
    b,
    c,
    pbc: [true, true, true],
    coordinateMode,
  };

  // Parse atoms
  const atoms: OpenQCAtom[] = [];
  let lineIdx = coordLine;

  for (let i = 0; i < elementNames.length; i++) {
    const element = elementNames[i];
    const count = atomCounts[i];

    for (let j = 0; j < count; j++) {
      if (lineIdx >= lines.length) {
        break;
      }
      const parts = lines[lineIdx].trim().split(/\s+/);
      lineIdx++;

      if (parts.length < 3) {
        continue;
      }

      const x = parseFloat(parts[0]);
      const y = parseFloat(parts[1]);
      const z = parseFloat(parts[2]);

      const atom: OpenQCAtom = { element, x, y, z };

      // Selective dynamics flags
      if (isSelective && parts.length >= 6) {
        const parseFlag = (s: string): boolean => s.toUpperCase().startsWith('T');
        atom.selectiveDynamics = [parseFlag(parts[3]), parseFlag(parts[4]), parseFlag(parts[5])];
      }

      atoms.push(atom);
    }
  }

  const structure: OpenQCStructure = {
    schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
    kind: 'periodic',
    name: filename ?? (comment || 'POSCAR'),
    atoms,
    cell,
    metadata: {
      source: {
        filename,
        software: 'vasp',
        parser: 'poscar',
      },
      provenance: {
        createdAt: new Date().toISOString(),
        warnings: [],
      },
    },
  };

  return structure;
}

// ---------------------------------------------------------------------------
// Convenience: create a minimal valid OpenQCStructure
// ---------------------------------------------------------------------------

/**
 * Create a minimal valid OpenQCStructure from a list of atoms.
 *
 * @param atoms - Array of atom data.
 * @param options - Optional name, kind, and cell.
 * @returns A valid OpenQCStructure.
 */
export function createOpenQCStructure(
  atoms: Array<{ element: string; x: number; y: number; z: number }>,
  options?: {
    name?: string;
    kind?: OpenQCStructureKind;
    cell?: OpenQCCell;
    charge?: number;
    multiplicity?: number;
  }
): OpenQCStructure {
  const kind = options?.kind ?? (options?.cell ? 'periodic' : 'molecule');
  const structure: OpenQCStructure = {
    schemaVersion: OPENQC_STRUCTURE_SCHEMA_VERSION,
    kind,
    name: options?.name,
    atoms: atoms.map(a => ({ element: a.element, x: a.x, y: a.y, z: a.z })),
    metadata: {
      charge: options?.charge,
      multiplicity: options?.multiplicity,
      provenance: {
        createdAt: new Date().toISOString(),
        warnings: [],
      },
    },
  };

  if (options?.cell) {
    structure.cell = options.cell;
  }

  return structure;
}
