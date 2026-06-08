/**
 * Runtime validation for OpenQCStructure DTOs.
 *
 * Provides a `validateOpenQCStructure` function that returns structured
 * diagnostics so callers can fail fast with clear error messages.
 *
 * @module structures/validation
 */

import type {
  OpenQCStructure,
  OpenQCAtom,
  OpenQCCell,
  OpenQCStructureKind,
  OpenQCStructureSchemaVersion,
} from './OpenQCStructure';
import { OPENQC_STRUCTURE_SCHEMA_VERSION } from './OpenQCStructure';

// ---------------------------------------------------------------------------
// Validation result
// ---------------------------------------------------------------------------

export interface ValidationResult {
  valid: boolean;
  errors: string[];
  warnings: string[];
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

const VALID_KINDS: ReadonlySet<string> = new Set<string>([
  'molecule',
  'periodic',
  'surface',
  'trajectory',
  'volumetric',
]);

const VALID_SCHEMA_VERSIONS: ReadonlySet<string> = new Set<string>([
  OPENQC_STRUCTURE_SCHEMA_VERSION,
]);

function isValidNumber(value: unknown): value is number {
  return typeof value === 'number' && isFinite(value) && !isNaN(value);
}

function isFiniteTuple3(value: unknown): value is [number, number, number] {
  return (
    Array.isArray(value) &&
    value.length === 3 &&
    isValidNumber(value[0]) &&
    isValidNumber(value[1]) &&
    isValidNumber(value[2])
  );
}

function validateAtom(atom: unknown, index: number, errors: string[], warnings: string[]): void {
  if (typeof atom !== 'object' || atom === null) {
    errors.push(`Atom ${index}: must be an object`);
    return;
  }

  const a = atom as Record<string, unknown>;

  // element
  if (typeof a.element !== 'string' || a.element.trim().length === 0) {
    errors.push(`Atom ${index}: missing or empty element symbol`);
  }

  // coordinates
  if (!isValidNumber(a.x)) {
    errors.push(`Atom ${index}: x coordinate is missing or not a finite number`);
  }
  if (!isValidNumber(a.y)) {
    errors.push(`Atom ${index}: y coordinate is missing or not a finite number`);
  }
  if (!isValidNumber(a.z)) {
    errors.push(`Atom ${index}: z coordinate is missing or not a finite number`);
  }

  // optional force tuple vs individual components
  if (a.force !== undefined && !isFiniteTuple3(a.force)) {
    errors.push(`Atom ${index}: force must be a 3-tuple of finite numbers`);
  }
  if (a.selectiveDynamics !== undefined) {
    if (
      !Array.isArray(a.selectiveDynamics) ||
      a.selectiveDynamics.length !== 3 ||
      a.selectiveDynamics.some(v => typeof v !== 'boolean')
    ) {
      errors.push(`Atom ${index}: selectiveDynamics must be a 3-tuple of booleans`);
    }
  }
}

function validateCell(cell: unknown, errors: string[]): void {
  if (typeof cell !== 'object' || cell === null) {
    errors.push('cell must be an object');
    return;
  }

  const c = cell as Record<string, unknown>;

  for (const key of ['a', 'b', 'c'] as const) {
    if (!isFiniteTuple3(c[key])) {
      errors.push(`cell.${key} must be a 3-tuple of finite numbers`);
    }
  }

  if (
    !Array.isArray(c.pbc) ||
    c.pbc.length !== 3 ||
    (c.pbc as unknown[]).some(v => typeof v !== 'boolean')
  ) {
    errors.push('cell.pbc must be a 3-tuple of booleans');
  }

  if (c.coordinateMode !== undefined) {
    if (c.coordinateMode !== 'cartesian' && c.coordinateMode !== 'fractional') {
      errors.push('cell.coordinateMode must be "cartesian" or "fractional"');
    }
  }

  // Warn about zero-length lattice vectors (degenerate cell)
  for (const key of ['a', 'b', 'c'] as const) {
    if (isFiniteTuple3(c[key])) {
      const v = c[key] as [number, number, number];
      const len = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
      if (len === 0) {
        errors.push(`cell.${key} has zero length -- degenerate lattice vector`);
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/**
 * Validate an OpenQCStructure at runtime.
 *
 * Returns a {@link ValidationResult} with `valid`, `errors`, and `warnings`.
 * A structure is considered valid when `errors` is empty.
 *
 * @param structure - The DTO to validate.
 * @param options   - Optional flags.
 * @param options.checkOverlap - Warn when two atoms are closer than tolerance (default true).
 * @param options.overlapTolerance - Overlap distance threshold in Angstroms (default 0.01).
 */
export function validateOpenQCStructure(
  structure: unknown,
  options?: { checkOverlap?: boolean; overlapTolerance?: number }
): ValidationResult {
  const errors: string[] = [];
  const warnings: string[] = [];
  const { checkOverlap = true, overlapTolerance = 0.01 } = options ?? {};

  if (typeof structure !== 'object' || structure === null) {
    return { valid: false, errors: ['structure must be a non-null object'], warnings };
  }

  const s = structure as Record<string, unknown>;

  // schemaVersion
  if (typeof s.schemaVersion !== 'string') {
    errors.push('schemaVersion is required and must be a string');
  } else if (!VALID_SCHEMA_VERSIONS.has(s.schemaVersion)) {
    errors.push(
      `Unsupported schemaVersion "${s.schemaVersion}". Supported: ${[...VALID_SCHEMA_VERSIONS].join(', ')}`
    );
  }

  // kind
  if (typeof s.kind !== 'string') {
    errors.push('kind is required and must be a string');
  } else if (!VALID_KINDS.has(s.kind)) {
    errors.push(
      `Invalid kind "${s.kind}". Must be one of: ${[...VALID_KINDS].join(', ')}`
    );
  }

  // atoms
  if (!Array.isArray(s.atoms)) {
    errors.push('atoms is required and must be an array');
  } else if (s.atoms.length === 0) {
    errors.push('atoms array must not be empty');
  } else {
    (s.atoms as unknown[]).forEach((atom, i) => validateAtom(atom, i, errors, warnings));
  }

  // cell (optional but must be valid when present)
  if (s.cell !== undefined) {
    validateCell(s.cell, errors);
  }

  // bonds (optional)
  if (s.bonds !== undefined) {
    if (!Array.isArray(s.bonds)) {
      errors.push('bonds must be an array when provided');
    } else {
      (s.bonds as unknown[]).forEach((bond, i) => {
        if (typeof bond !== 'object' || bond === null) {
          errors.push(`Bond ${i}: must be an object`);
          return;
        }
        const b = bond as Record<string, unknown>;
        if (typeof b.from !== 'number' || typeof b.to !== 'number') {
          errors.push(`Bond ${i}: "from" and "to" must be numbers`);
        }
        if (b.order !== undefined && typeof b.order !== 'number') {
          errors.push(`Bond ${i}: "order" must be a number when provided`);
        }
      });
    }
  }

  // frames (optional)
  if (s.frames !== undefined) {
    if (!Array.isArray(s.frames)) {
      errors.push('frames must be an array when provided');
    }
    // Recursively validate each frame (lightweight -- one level deep to avoid
    // deep recursion on malformed data).
    (s.frames as unknown[]).forEach((frame, i) => {
      if (typeof frame !== 'object' || frame === null) {
        errors.push(`Frame ${i}: must be an object`);
      } else {
        const f = frame as Record<string, unknown>;
        if (!Array.isArray(f.atoms) || (f.atoms as unknown[]).length === 0) {
          warnings.push(`Frame ${i}: has no atoms`);
        }
      }
    });
  }

  // metadata (optional)
  if (s.metadata !== undefined) {
    if (typeof s.metadata !== 'object' || s.metadata === null) {
      errors.push('metadata must be an object when provided');
    }
  }

  // Duplicate-atom proximity check (warning only)
  if (checkOverlap && Array.isArray(s.atoms) && s.atoms.length > 1) {
    const atoms = s.atoms as OpenQCAtom[];
    for (let i = 0; i < atoms.length; i++) {
      for (let j = i + 1; j < atoms.length; j++) {
        const dx = atoms[i].x - atoms[j].x;
        const dy = atoms[i].y - atoms[j].y;
        const dz = atoms[i].z - atoms[j].z;
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (dist < overlapTolerance) {
          warnings.push(`Atoms ${i} and ${j} are very close (${dist.toFixed(4)} A)`);
        }
      }
    }
  }

  return { valid: errors.length === 0, errors, warnings };
}
