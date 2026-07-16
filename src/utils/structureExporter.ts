/**
 * Structure Exporter - Export modified structures to various formats
 *
 * Phase 2: Week 7-8 - Export modified structures back to input format
 */

import * as vscode from 'vscode';
import { Atom } from '../visualizers/types';
import type { OpenQCBond, OpenQCCell } from '../structures/OpenQCStructure';

export type ExportFormat =
  | 'vasp'
  | 'cp2k'
  | 'qe'
  | 'gaussian'
  | 'orca'
  | 'nwchem'
  | 'gamess'
  | 'lammps'
  | 'xyz'
  | 'extxyz'
  | 'pdb'
  | 'cif';

export interface ExportOptions {
  format: ExportFormat;
  outputPath?: string;
  preserveComments?: boolean;
  preserveMetadata?: boolean;
  overwrite?: boolean;
}

export interface StructureData {
  atoms: Atom[];
  bonds?: OpenQCBond[];
  cell?: LatticeParameters | OpenQCCell;
  pbc?: [boolean, boolean, boolean];
  metadata?: Record<string, any>;
}

export interface LatticeParameters {
  a: number;
  b: number;
  c: number;
  alpha: number;
  beta: number;
  gamma: number;
}

export interface ExportResult {
  success: boolean;
  outputPath?: string;
  content?: string;
  error?: string;
  warnings?: string[];
}

export class StructureExporter {
  constructor(_context?: vscode.ExtensionContext) {}

  public async exportStructure(
    structure: StructureData,
    options: ExportOptions
  ): Promise<ExportResult> {
    if (!structure.atoms || structure.atoms.length === 0) {
      return { success: false, error: 'No atoms to export' };
    }

    const contentResult = this.renderStructure(structure, options.format);
    if (!contentResult.success || !contentResult.content) {
      return contentResult;
    }

    let outputPath = options.outputPath;
    if (!outputPath) {
      const uri = await vscode.window.showSaveDialog({
        defaultUri: vscode.Uri.file(`structure.${getFormatExtensions(options.format)[0]}`),
        saveLabel: 'Export Structure',
        filters: {
          [StructureExporter.getFormatDisplayName(options.format)]: getFormatExtensions(
            options.format
          ),
        },
      });
      if (!uri) {
        return { success: false, error: 'Export cancelled' };
      }
      outputPath = uri.fsPath;
    }

    await vscode.workspace.fs.writeFile(
      vscode.Uri.file(outputPath),
      Buffer.from(contentResult.content)
    );
    vscode.window.showInformationMessage(`Structure exported to ${outputPath}`);

    return {
      success: true,
      outputPath,
      content: contentResult.content,
      warnings: contentResult.warnings,
    };
  }

  public async overwriteStructureFile(
    structure: StructureData,
    outputPath: string
  ): Promise<ExportResult> {
    const format = StructureExporter.inferNativeFormatFromPath(outputPath);
    if (!format) {
      return {
        success: false,
        error:
          'Cannot infer a native writable structure format for this source file. Export to XYZ, Extended XYZ, PDB, POSCAR, or CIF instead.',
      };
    }

    return this.exportStructure(structure, {
      format,
      outputPath,
      overwrite: true,
    });
  }

  private renderStructure(structure: StructureData, format: ExportFormat): ExportResult {
    switch (format) {
      case 'xyz':
        return {
          success: true,
          content: renderXyz(structure),
          warnings: bondTopologyWarnings('XYZ', structure),
        };
      case 'extxyz':
        return {
          success: true,
          content: renderExtxyz(structure),
          warnings: bondTopologyWarnings('Extended XYZ', structure),
        };
      case 'pdb':
        return renderPdb(structure);
      case 'vasp':
        return renderPoscar(structure);
      case 'cif':
        return renderCif(structure);
      default:
        return {
          success: false,
          error: `${StructureExporter.getFormatDisplayName(format)} export requires the ASE backend and is not implemented by the native exporter yet.`,
        };
    }
  }

  public static getFormatDisplayName(format: ExportFormat): string {
    const nameMap: Record<ExportFormat, string> = {
      vasp: 'VASP POSCAR',
      cp2k: 'CP2K Input',
      qe: 'Quantum ESPRESSO',
      gaussian: 'Gaussian Input',
      orca: 'ORCA Input',
      nwchem: 'NWChem Input',
      gamess: 'GAMESS Input',
      lammps: 'LAMMPS Data',
      xyz: 'XYZ',
      extxyz: 'Extended XYZ',
      pdb: 'PDB',
      cif: 'CIF',
    };
    return nameMap[format];
  }

  public static getSupportedFormats(): ExportFormat[] {
    return [
      'vasp',
      'cp2k',
      'qe',
      'gaussian',
      'orca',
      'nwchem',
      'gamess',
      'lammps',
      'xyz',
      'extxyz',
      'pdb',
      'cif',
    ];
  }

  public static getNativeFormats(): ExportFormat[] {
    return ['xyz', 'extxyz', 'pdb', 'vasp', 'cif'];
  }

  public static inferNativeFormatFromPath(filePath: string): ExportFormat | undefined {
    const normalized = filePath.replace(/\\/g, '/');
    const basename = normalized.split('/').pop() ?? normalized;
    const upperBasename = basename.toUpperCase();
    const lowerBasename = basename.toLowerCase();

    if (upperBasename === 'POSCAR' || upperBasename === 'CONTCAR') {
      return 'vasp';
    }
    if (lowerBasename.endsWith('.extxyz')) {
      return 'extxyz';
    }
    if (lowerBasename.endsWith('.xyz')) {
      return 'xyz';
    }
    if (lowerBasename.endsWith('.pdb')) {
      return 'pdb';
    }
    if (lowerBasename.endsWith('.cif')) {
      return 'cif';
    }

    return undefined;
  }
}

function bondTopologyWarnings(formatName: string, structure: StructureData): string[] {
  if (!Array.isArray(structure.bonds) || structure.bonds.length === 0) {
    return [];
  }
  return [
    `${formatName} export does not preserve edited bond topology or bond order; use PDB for basic CONECT connectivity.`,
  ];
}

function renderXyz(structure: StructureData): string {
  const atoms = cartesianAtoms(structure);
  return (
    [
      String(structure.atoms.length),
      'OpenQC export',
      ...atoms.map(
        atom =>
          `${atom.element} ${formatCoord(atom.x)} ${formatCoord(atom.y)} ${formatCoord(atom.z)}`
      ),
    ].join('\n') + '\n'
  );
}

function renderExtxyz(structure: StructureData): string {
  const vectors = getCellVectors(structure.cell);
  const atoms = cartesianAtoms(structure);
  const lattice = vectors
    ? ` Lattice="${vectors.flat().map(formatCoord).join(' ')}" Properties=species:S:1:pos:R:3`
    : ' Properties=species:S:1:pos:R:3';
  return (
    [
      String(structure.atoms.length),
      `OpenQC export${lattice}`,
      ...atoms.map(
        atom =>
          `${atom.element} ${formatCoord(atom.x)} ${formatCoord(atom.y)} ${formatCoord(atom.z)}`
      ),
    ].join('\n') + '\n'
  );
}

function renderPdb(structure: StructureData): ExportResult {
  const lines = cartesianAtoms(structure).map((atom, index) => {
    const serial = String(index + 1).padStart(5);
    const element = atom.element.padStart(2);
    return `HETATM${serial} ${element.padEnd(4)} MOL     1    ${formatPdbCoord(atom.x)}${formatPdbCoord(atom.y)}${formatPdbCoord(atom.z)}  1.00  0.00          ${element}`;
  });
  const conect = renderPdbConect(structure.bonds, structure.atoms.length);
  return {
    success: true,
    content: [...lines, ...conect.lines, 'END'].join('\n') + '\n',
    warnings: conect.warnings,
  };
}

function renderPdbConect(
  bonds: OpenQCBond[] | undefined,
  atomCount: number
): { lines: string[]; warnings: string[] } {
  if (!Array.isArray(bonds) || bonds.length === 0) {
    return { lines: [], warnings: [] };
  }

  const lines: string[] = [];
  const warnings: string[] = [];
  const seen = new Set<string>();
  for (const [index, bond] of bonds.entries()) {
    if (!isValidBondIndex(bond.from, atomCount) || !isValidBondIndex(bond.to, atomCount)) {
      warnings.push(
        `Skipped PDB CONECT bond ${index}: atom indices ${bond.from}-${bond.to} are out of range for ${atomCount} atoms`
      );
      continue;
    }
    if (bond.from === bond.to) {
      warnings.push(`Skipped PDB CONECT bond ${index}: self bonds cannot be represented`);
      continue;
    }
    const from = Math.min(bond.from, bond.to);
    const to = Math.max(bond.from, bond.to);
    const key = `${from}:${to}`;
    if (seen.has(key)) {
      continue;
    }
    seen.add(key);
    lines.push(`CONECT${String(from + 1).padStart(5)}${String(to + 1).padStart(5)}`);
  }
  return { lines, warnings };
}

function isValidBondIndex(index: number, atomCount: number): boolean {
  return Number.isInteger(index) && index >= 0 && index < atomCount;
}

function renderPoscar(structure: StructureData): ExportResult {
  const vectors = getCellVectors(structure.cell);
  if (!vectors) {
    return {
      success: false,
      error: 'VASP POSCAR export requires unit-cell vectors or lattice parameters',
    };
  }

  const groups = groupAtomsByElement(structure.atoms);
  const coordinateMode = getCoordinateMode(structure.cell);
  const hasSelectiveDynamics = structure.atoms.some(atom => hasSelectiveDynamicsFlags(atom));
  const lines = [
    'OpenQC export',
    '1.0',
    ...vectors.map(vector => vector.map(formatCoord).join(' ')),
    groups.map(group => group.element).join(' '),
    groups.map(group => String(group.atoms.length)).join(' '),
  ];

  if (hasSelectiveDynamics) {
    lines.push('Selective dynamics');
  }
  lines.push(coordinateMode === 'fractional' ? 'Direct' : 'Cartesian');

  for (const group of groups) {
    for (const atom of group.atoms) {
      const coords = `${formatCoord(atom.x)} ${formatCoord(atom.y)} ${formatCoord(atom.z)}`;
      lines.push(hasSelectiveDynamics ? `${coords} ${formatSelectiveDynamics(atom)}` : coords);
    }
  }

  return {
    success: true,
    content: lines.join('\n') + '\n',
    warnings: bondTopologyWarnings('VASP POSCAR', structure),
  };
}

function hasSelectiveDynamicsFlags(atom: Atom): boolean {
  return (
    Array.isArray(atom.selectiveDynamics) &&
    atom.selectiveDynamics.length === 3 &&
    atom.selectiveDynamics.every(value => typeof value === 'boolean')
  );
}

function formatSelectiveDynamics(atom: Atom): string {
  const flags: [boolean, boolean, boolean] = hasSelectiveDynamicsFlags(atom)
    ? (atom.selectiveDynamics as [boolean, boolean, boolean])
    : [true, true, true];
  return flags.map(value => (value ? 'T' : 'F')).join(' ');
}

function renderCif(structure: StructureData): ExportResult {
  const vectors = getCellVectors(structure.cell);
  if (!vectors) {
    return { success: false, error: 'CIF export requires unit-cell vectors or lattice parameters' };
  }

  const params = latticeParametersFromVectors(vectors);
  const atoms = cartesianAtoms(structure);
  const lines = [
    'data_openqc_export',
    `_cell_length_a ${formatCoord(params.a)}`,
    `_cell_length_b ${formatCoord(params.b)}`,
    `_cell_length_c ${formatCoord(params.c)}`,
    `_cell_angle_alpha ${formatCoord(params.alpha)}`,
    `_cell_angle_beta ${formatCoord(params.beta)}`,
    `_cell_angle_gamma ${formatCoord(params.gamma)}`,
    'loop_',
    '_atom_site_label',
    '_atom_site_type_symbol',
    '_atom_site_Cartn_x',
    '_atom_site_Cartn_y',
    '_atom_site_Cartn_z',
  ];

  atoms.forEach((atom, index) => {
    const label = `${atom.element}${index + 1}`;
    lines.push(
      `${label} ${atom.element} ${formatCoord(atom.x)} ${formatCoord(atom.y)} ${formatCoord(atom.z)}`
    );
  });

  return {
    success: true,
    content: lines.join('\n') + '\n',
    warnings: bondTopologyWarnings('CIF', structure),
  };
}

function cartesianAtoms(structure: StructureData): Atom[] {
  const vectors = getCellVectors(structure.cell);
  if (!vectors || getCoordinateMode(structure.cell) !== 'fractional') {
    return structure.atoms;
  }

  return structure.atoms.map(atom => {
    const [x, y, z] = fractionalToCartesian([atom.x, atom.y, atom.z], vectors);
    return { ...atom, x, y, z };
  });
}

function getCoordinateMode(cell: StructureData['cell']): 'cartesian' | 'fractional' {
  return (cell as OpenQCCell | undefined)?.coordinateMode === 'fractional'
    ? 'fractional'
    : 'cartesian';
}

function fractionalToCartesian(
  fractional: [number, number, number],
  vectors: [[number, number, number], [number, number, number], [number, number, number]]
): [number, number, number] {
  const [u, v, w] = fractional;
  const [a, b, c] = vectors;
  return [
    u * a[0] + v * b[0] + w * c[0],
    u * a[1] + v * b[1] + w * c[1],
    u * a[2] + v * b[2] + w * c[2],
  ];
}

function getCellVectors(
  cell: StructureData['cell']
): [[number, number, number], [number, number, number], [number, number, number]] | undefined {
  if (!cell) {
    return undefined;
  }

  if (Array.isArray((cell as OpenQCCell).a)) {
    const openqcCell = cell as OpenQCCell;
    return [openqcCell.a, openqcCell.b, openqcCell.c];
  }

  const params = cell as LatticeParameters;
  const alpha = degreesToRadians(params.alpha);
  const beta = degreesToRadians(params.beta);
  const gamma = degreesToRadians(params.gamma);
  const a: [number, number, number] = [params.a, 0, 0];
  const b: [number, number, number] = [params.b * Math.cos(gamma), params.b * Math.sin(gamma), 0];
  const cx = params.c * Math.cos(beta);
  const cy = (params.c * (Math.cos(alpha) - Math.cos(beta) * Math.cos(gamma))) / Math.sin(gamma);
  const cz = Math.sqrt(Math.max(params.c ** 2 - cx ** 2 - cy ** 2, 0));
  return [a, b, [cx, cy, cz]];
}

function latticeParametersFromVectors(
  vectors: [[number, number, number], [number, number, number], [number, number, number]]
): LatticeParameters {
  const [a, b, c] = vectors;
  return {
    a: norm(a),
    b: norm(b),
    c: norm(c),
    alpha: angleDegrees(b, c),
    beta: angleDegrees(a, c),
    gamma: angleDegrees(a, b),
  };
}

function groupAtomsByElement(atoms: Atom[]): Array<{ element: string; atoms: Atom[] }> {
  const groups = new Map<string, Atom[]>();
  for (const atom of atoms) {
    const group = groups.get(atom.element) ?? [];
    group.push(atom);
    groups.set(atom.element, group);
  }
  return Array.from(groups, ([element, atoms]) => ({ element, atoms }));
}

function norm(vector: [number, number, number]): number {
  return Math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2);
}

function angleDegrees(a: [number, number, number], b: [number, number, number]): number {
  const dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const cosine = Math.min(1, Math.max(-1, dot / (norm(a) * norm(b))));
  return (Math.acos(cosine) * 180) / Math.PI;
}

function degreesToRadians(value: number): number {
  return (value * Math.PI) / 180;
}

function formatCoord(value: number): string {
  return Number.isFinite(value) ? value.toFixed(6) : '0.000000';
}

function formatPdbCoord(value: number): string {
  return value.toFixed(3).padStart(8);
}

function getFormatExtensions(format: ExportFormat): string[] {
  const extensionMap: Record<ExportFormat, string[]> = {
    vasp: ['POSCAR', 'CONTCAR'],
    cp2k: ['inp'],
    qe: ['in', 'pw.in'],
    gaussian: ['com', 'gjf'],
    orca: ['inp'],
    nwchem: ['nw', 'nwinp'],
    gamess: ['inp'],
    lammps: ['lmp', 'data'],
    xyz: ['xyz'],
    extxyz: ['extxyz'],
    pdb: ['pdb'],
    cif: ['cif'],
  };
  return extensionMap[format];
}
