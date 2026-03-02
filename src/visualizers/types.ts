/**
 * Type definitions for Three.js molecular visualization
 *
 * Phase 2: 3D Visualization - Types and Interfaces
 */

/**
 * Represents an atom with its 3D coordinates and element type
 */
export interface Atom {
  element: string;
  x: number;
  y: number;
  z: number;
}

/**
 * Represents a bond between two atoms
 */
export interface Bond {
  atomIndex1: number;
  atomIndex2: number;
  length: number;
  order?: number; // single, double, triple bond
}

/**
 * Lattice vectors for periodic systems
 */
export interface LatticeVectors {
  a: [number, number, number];
  b: [number, number, number];
  c: [number, number, number];
}

/**
 * Molecular structure data
 */
export interface MolecularStructure {
  atoms: Atom[];
  bonds?: Bond[];
  lattice?: LatticeVectors;
  name?: string;
}

/**
 * Color definition for elements
 */
export interface ElementColor {
  element: string;
  color: string; // Hex color code
  rgb: { r: number; g: number; b: number };
}

/**
 * Atomic radii for visualization
 */
export interface AtomicRadii {
  covalent: number;
  vanDerWaals: number;
}

/**
 * Representation modes for molecular visualization
 */
export type RepresentationMode =
  | 'ball-and-stick'
  | 'space-filling'
  | 'wireframe'
  | 'stick';

/**
 * Camera position and orientation
 */
export interface CameraState {
  position: { x: number; y: number; z: number };
  target: { x: number; y: number; z: number };
  zoom: number;
}

/**
 * Render result information
 */
export interface RenderResult {
  success: boolean;
  atomCount: number;
  bondCount?: number;
  hasCell?: boolean;
  renderTime?: number;
}

/**
 * Export options for snapshots
 */
export interface ExportOptions {
  format: 'png' | 'jpeg' | 'webp';
  quality?: number; // 0-1 for jpeg/webp
  width?: number;
  height?: number;
  backgroundColor?: string;
}

/**
 * Visualization configuration
 */
export interface VisualizationConfig {
  representationMode: RepresentationMode;
  showBonds: boolean;
  showUnitCell: boolean;
  showAxes: boolean;
  atomScale: number;
  bondRadius: number;
  bgColor: string;
  antialiasing: boolean;
}

/**
 * Three.js renderer state
 */
export interface RendererState {
  atoms: Atom[];
  bonds: Bond[];
  camera: CameraState;
  config: VisualizationConfig;
}

/**
 * Element colors following CPK coloring convention
 */
export const ELEMENT_COLORS: Record<string, string> = {
  H: '#FFFFFF',
  C: '#909090',
  N: '#3050F8',
  O: '#FF0D0D',
  F: '#90E050',
  P: '#FF8000',
  S: '#FFFF30',
  Cl: '#1FF01F',
  Br: '#A62929',
  I: '#940094',
  // Metals
  Li: '#CC80FF',
  Na: '#AB5CF2',
  K: '#8F40D4',
  Mg: '#8AFF00',
  Ca: '#3DFF00',
  Fe: '#E06633',
  Cu: '#C88033',
  Zn: '#7D80B0',
  Ag: '#C0C0C0',
  Au: '#FFD123',
  Hg: '#B8B8D0',
  Pb: '#575961',
  // Default for unknown elements
  default: '#FF1493',
};

/**
 * Covalent radii in Angstroms (for bond detection)
 */
export const COVALENT_RADII: Record<string, number> = {
  H: 0.31,
  C: 0.76,
  N: 0.71,
  O: 0.66,
  F: 0.57,
  P: 1.07,
  S: 1.05,
  Cl: 1.02,
  Br: 1.2,
  I: 1.39,
  Li: 1.28,
  Na: 1.66,
  K: 2.03,
  Mg: 1.41,
  Ca: 1.76,
  Fe: 1.32,
  Cu: 1.32,
  Zn: 1.22,
  Ag: 1.45,
  Au: 1.36,
  Hg: 1.32,
  Pb: 1.46,
  default: 1.5,
};

/**
 * Van der Waals radii in Angstroms (for space-filling mode)
 */
export const VDW_RADII: Record<string, number> = {
  H: 1.2,
  C: 1.7,
  N: 1.55,
  O: 1.52,
  F: 1.47,
  P: 1.8,
  S: 1.8,
  Cl: 1.75,
  Br: 1.85,
  I: 1.98,
  Li: 1.82,
  Na: 2.27,
  K: 2.75,
  Mg: 1.73,
  Ca: 2.31,
  Fe: 1.94,
  Cu: 1.87,
  Zn: 2.1,
  Ag: 2.1,
  Au: 2.14,
  Hg: 2.23,
  Pb: 2.32,
  default: 2.0,
};
