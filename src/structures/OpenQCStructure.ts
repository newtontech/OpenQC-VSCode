/**
 * OpenQCStructure - Single source-of-truth DTO for molecular/crystal structures.
 *
 * Every parser, viewer, converter, and agent tool in OpenQC should use this
 * canonical representation instead of ad-hoc atom arrays or XYZ strings.
 *
 * @module structures/OpenQCStructure
 */

// ---------------------------------------------------------------------------
// Schema version constant
// ---------------------------------------------------------------------------

export const OPENQC_STRUCTURE_SCHEMA_VERSION = 'openqc.structure.v1' as const;
export type OpenQCStructureSchemaVersion = typeof OPENQC_STRUCTURE_SCHEMA_VERSION;

// ---------------------------------------------------------------------------
// Structure kind
// ---------------------------------------------------------------------------

export type OpenQCStructureKind = 'molecule' | 'periodic' | 'surface' | 'trajectory' | 'volumetric';

// ---------------------------------------------------------------------------
// Atom
// ---------------------------------------------------------------------------

export interface OpenQCAtom {
  /** Element symbol (e.g. "Fe", "C"). */
  element: string;
  /** Cartesian x coordinate in Angstroms. */
  x: number;
  /** Cartesian y coordinate in Angstroms. */
  y: number;
  /** Cartesian z coordinate in Angstroms. */
  z: number;
  /** Cartesian force x component (optional, eV/Angstrom). */
  fx?: number;
  /** Cartesian force y component (optional, eV/Angstrom). */
  fy?: number;
  /** Cartesian force z component (optional, eV/Angstrom). */
  fz?: number;
  /** Atomic charge (optional). */
  charge?: number;
  /** Magnetic moment (optional, mu_B). */
  magmom?: number;
  /** Force as a 3-tuple (optional, alternative to fx/fy/fz). */
  force?: [number, number, number];
  /** Per-atom label (optional, e.g. "C1"). */
  label?: string;
  /** Selective dynamics flags (VASP-style). */
  selectiveDynamics?: [boolean, boolean, boolean];
}

// ---------------------------------------------------------------------------
// Cell / lattice
// ---------------------------------------------------------------------------

export interface OpenQCCell {
  /** Lattice vector a in Angstroms. */
  a: [number, number, number];
  /** Lattice vector b in Angstroms. */
  b: [number, number, number];
  /** Lattice vector c in Angstroms. */
  c: [number, number, number];
  /** Periodic boundary condition flags. */
  pbc: [boolean, boolean, boolean];
  /** Coordinate mode for atom positions. */
  coordinateMode?: 'cartesian' | 'fractional';
}

// ---------------------------------------------------------------------------
// Bond
// ---------------------------------------------------------------------------

export interface OpenQCBond {
  /** Zero-based index of the first atom. */
  from: number;
  /** Zero-based index of the second atom. */
  to: number;
  /** Bond order (1=single, 2=double, 3=triple). */
  order?: number;
  /** Periodic image offset for bonds across cell boundaries. */
  periodicImage?: [number, number, number];
}

// ---------------------------------------------------------------------------
// Structure metadata
// ---------------------------------------------------------------------------

export interface OpenQCStructureMetadata {
  /** Charge of the system. */
  charge?: number;
  /** Spin multiplicity. */
  multiplicity?: number;
  /** Total energy (optional, eV). */
  energyEv?: number;
  /** Source file information. */
  source?: {
    uri?: string;
    filename?: string;
    software?: string;
    parser?: string;
  };
  /** Provenance information. */
  provenance?: {
    createdAt: string;
    openqcVersion?: string;
    warnings: string[];
  };
}

// ---------------------------------------------------------------------------
// OpenQCStructure (the top-level DTO)
// ---------------------------------------------------------------------------

export interface OpenQCStructure {
  /** Schema version for forward-compatible migration. */
  schemaVersion: OpenQCStructureSchemaVersion;
  /** Discriminator indicating the type of structure. */
  kind: OpenQCStructureKind;
  /** Human-readable name for the structure. */
  name?: string;
  /** Atom list -- always Cartesian unless cell.coordinateMode says otherwise. */
  atoms: OpenQCAtom[];
  /** Unit cell for periodic systems. Omit for isolated molecules. */
  cell?: OpenQCCell;
  /** Optional bond list. */
  bonds?: OpenQCBond[];
  /** Frames for trajectory support (each frame is a full OpenQCStructure). */
  frames?: OpenQCStructure[];
  /** Arbitrary metadata. */
  metadata?: OpenQCStructureMetadata;
}
