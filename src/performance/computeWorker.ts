/**
 * Compute Worker - WebWorker for Heavy Computations
 *
 * Offloads computationally intensive operations to background threads
 * to prevent UI blocking in VSCode.
 */

import { ASEAtoms } from '../ase/ASEConverter';
import {
  convertParameterValue,
  getParameterMappings,
  ParameterMapping,
} from '../utils/migration/params';

// Worker message types
export enum WorkerMessageType {
  PARSE_STRUCTURE = 'PARSE_STRUCTURE',
  CONVERT_FORMAT = 'CONVERT_FORMAT',
  VALIDATE_STRUCTURE = 'VALIDATE_STRUCTURE',
  CALCULATE_PROPERTIES = 'CALCULATE_PROPERTIES',
  MIGRATE_PARAMETERS = 'MIGRATE_PARAMETERS',
}

export interface WorkerMessage {
  type: WorkerMessageType;
  id: string;
  payload: any;
}

export interface WorkerResponse {
  id: string;
  success: boolean;
  result?: any;
  error?: string;
  duration: number;
}

/**
 * Structure parsing task
 */
interface ParseStructurePayload {
  content: string;
  format: string;
  options?: Record<string, any>;
}

/**
 * Format conversion task
 */
interface ConvertFormatPayload {
  atoms: ASEAtoms;
  targetFormat: string;
  options?: Record<string, any>;
}

/**
 * Structure validation task
 */
interface ValidateStructurePayload {
  atoms: ASEAtoms;
  checks: string[];
}

/**
 * Property calculation task
 */
interface CalculatePropertiesPayload {
  atoms: ASEAtoms;
  properties: string[];
}

/**
 * Parameter migration task
 */
interface MigrateParametersPayload {
  sourceFormat: string;
  targetFormat: string;
  parameters: Record<string, any>;
}

/**
 * Worker implementation
 */
class ComputeWorker {
  private startTime: number = 0;

  /**
   * Process incoming message
   */
  async processMessage(message: WorkerMessage): Promise<WorkerResponse> {
    this.startTime = Date.now();

    try {
      let result: any;

      switch (message.type) {
        case WorkerMessageType.PARSE_STRUCTURE:
          result = await this.parseStructure(message.payload as ParseStructurePayload);
          break;

        case WorkerMessageType.CONVERT_FORMAT:
          result = await this.convertFormat(message.payload as ConvertFormatPayload);
          break;

        case WorkerMessageType.VALIDATE_STRUCTURE:
          result = await this.validateStructure(message.payload as ValidateStructurePayload);
          break;

        case WorkerMessageType.CALCULATE_PROPERTIES:
          result = await this.calculateProperties(message.payload as CalculatePropertiesPayload);
          break;

        case WorkerMessageType.MIGRATE_PARAMETERS:
          result = await this.migrateParameters(message.payload as MigrateParametersPayload);
          break;

        default:
          throw new Error(`Unknown message type: ${message.type}`);
      }

      return {
        id: message.id,
        success: true,
        result,
        duration: Date.now() - this.startTime,
      };
    } catch (error) {
      return {
        id: message.id,
        success: false,
        error: error instanceof Error ? error.message : String(error),
        duration: Date.now() - this.startTime,
      };
    }
  }

  /**
   * Parse structure from text content
   */
  private async parseStructure(payload: ParseStructurePayload): Promise<ASEAtoms> {
    const { content, format } = payload;

    if (!content || content.trim().length === 0) {
      throw new Error('Empty content provided');
    }

    const normalizedFormat = String(format || '').toLowerCase();
    if (['xyz', 'extxyz'].includes(normalizedFormat)) {
      return this.parseXyz(content, normalizedFormat);
    }
    if (['vasp', 'poscar', 'contcar'].includes(normalizedFormat)) {
      return this.parsePoscar(content, normalizedFormat);
    }

    throw new Error(`Unsupported structure format for worker-local parsing: ${format}`);
  }

  /**
   * Convert atoms to target format
   */
  private async convertFormat(payload: ConvertFormatPayload): Promise<string> {
    const { atoms, targetFormat } = payload;
    const normalizedTarget = String(targetFormat || '').toLowerCase();

    this.validateAtoms(atoms);

    if (['xyz', 'extxyz'].includes(normalizedTarget)) {
      return this.formatXyz(atoms, normalizedTarget === 'extxyz');
    }
    if (['vasp', 'poscar', 'contcar'].includes(normalizedTarget)) {
      return this.formatPoscar(atoms);
    }

    throw new Error(`Unsupported worker conversion target format: ${targetFormat}`);
  }

  /**
   * Validate structure
   */
  private async validateStructure(payload: ValidateStructurePayload): Promise<Record<string, any>> {
    const { atoms, checks } = payload;
    const results: Record<string, any> = {};
    this.validateAtoms(atoms);

    for (const check of checks) {
      switch (check) {
        case 'bond_lengths':
          results[check] = this.checkBondLengths(atoms);
          break;

        case 'cell_consistency':
          results[check] = this.checkCellConsistency(atoms);
          break;

        case 'atom_overlap':
          results[check] = this.checkAtomOverlap(atoms);
          break;

        case 'charge_neutrality':
          results[check] = this.checkChargeNeutrality(atoms);
          break;

        default:
          results[check] = { valid: true, warnings: [`Unknown check: ${check}`] };
      }
    }

    return results;
  }

  /**
   * Calculate molecular properties
   */
  private async calculateProperties(
    payload: CalculatePropertiesPayload
  ): Promise<Record<string, any>> {
    const { atoms, properties } = payload;
    const results: Record<string, any> = {};
    this.validateAtoms(atoms);

    for (const prop of properties) {
      switch (prop) {
        case 'center_of_mass':
          results[prop] = this.calculateCenterOfMass(atoms);
          break;

        case 'moment_of_inertia':
          results[prop] = this.calculateMomentOfInertia(atoms);
          break;

        case 'bounding_box':
          results[prop] = this.calculateBoundingBox(atoms);
          break;

        case 'atom_count':
          results[prop] = atoms.chemical_symbols.length;
          break;

        default:
          results[prop] = null;
      }
    }

    return results;
  }

  /**
   * Migrate parameters between formats
   */
  private async migrateParameters(payload: MigrateParametersPayload): Promise<Record<string, any>> {
    const { sourceFormat, targetFormat, parameters } = payload;
    const source = String(sourceFormat || '').toLowerCase();
    const target = String(targetFormat || '').toLowerCase();
    const migrated: Record<string, any> = {};
    const warnings: string[] = [];
    const notes = [`Migrated from ${source} to ${target}`];
    const unmapped: string[] = [];

    if (!parameters || typeof parameters !== 'object' || Array.isArray(parameters)) {
      throw new Error('parameters must be a key-value object');
    }

    const mappings = getParameterMappings(source, target);
    if (mappings.length === 0) {
      warnings.push(`No worker-local parameter mappings for ${source} -> ${target}`);
      return {
        migrated,
        warnings,
        notes,
        unmapped: Object.keys(parameters),
        supported: false,
      };
    }

    const mappingLookup = new Map<string, ParameterMapping>();
    mappings.forEach(mapping => mappingLookup.set(mapping.sourceParam.toUpperCase(), mapping));

    Object.entries(parameters).forEach(([key, value]) => {
      const baseKey = this.getBaseParameterName(key);
      const mapping =
        mappingLookup.get(key.toUpperCase()) || mappingLookup.get(baseKey.toUpperCase());

      if (!mapping) {
        unmapped.push(key);
        return;
      }

      try {
        const convertedValue = convertParameterValue(mapping, value);
        if (
          mapping.unitFactor !== undefined &&
          typeof convertedValue === 'number' &&
          !Number.isFinite(convertedValue)
        ) {
          throw new Error(`non-numeric value cannot be converted with unit factor`);
        }
        migrated[mapping.targetParam] = this.cleanConvertedValue(convertedValue);
        if (mapping.unitFactor !== undefined && mapping.unitFactor !== 1) {
          notes.push(
            `${key} -> ${mapping.targetParam} converted with factor ${mapping.unitFactor}`
          );
        }
      } catch (error) {
        warnings.push(
          `Failed to map ${key} to ${mapping.targetParam}: ${
            error instanceof Error ? error.message : String(error)
          }`
        );
      }
    });

    if (unmapped.length > 0) {
      warnings.push(
        `${unmapped.length} parameters could not be mapped: ${unmapped
          .slice(0, 5)
          .join(', ')}${unmapped.length > 5 ? '...' : ''}`
      );
    }

    return {
      migrated,
      warnings,
      notes,
      unmapped,
      supported: true,
    };
  }

  /**
   * Check bond lengths
   */
  private checkBondLengths(atoms: ASEAtoms): any {
    const warnings: string[] = [];
    let minDist = Infinity;
    let maxDist = 0;

    // Simple distance check (would be more sophisticated in production)
    for (let i = 0; i < Math.min(atoms.chemical_symbols.length, 100); i++) {
      for (let j = i + 1; j < Math.min(atoms.chemical_symbols.length, 100); j++) {
        const dist = this.distance(atoms.positions[i], atoms.positions[j]);
        if (dist < 0.5) {
          warnings.push(`Atoms ${i}-${j} too close: ${dist.toFixed(3)} Å`);
        }
        minDist = Math.min(minDist, dist);
        maxDist = Math.max(maxDist, dist);
      }
    }

    return {
      valid: warnings.length === 0,
      warnings,
      minDistance: minDist,
      maxDistance: maxDist,
    };
  }

  /**
   * Check cell consistency
   */
  private checkCellConsistency(atoms: ASEAtoms): any {
    if (!atoms.pbc.some(p => p)) {
      return { valid: true, warnings: ['Non-periodic system'] };
    }
    if (!atoms.cell) {
      return { valid: false, warnings: ['Periodic system is missing unit-cell vectors'] };
    }

    const warnings: string[] = [];

    // Check cell vectors
    for (let i = 0; i < 3; i++) {
      const length = this.distance(atoms.cell![i], [0, 0, 0]);
      if (length < 1.0) {
        warnings.push(`Cell vector ${i + 1} too short: ${length.toFixed(3)} Å`);
      }
    }

    return {
      valid: warnings.length === 0,
      warnings,
    };
  }

  /**
   * Check for atom overlap
   */
  private checkAtomOverlap(atoms: ASEAtoms): any {
    const overlaps: Array<[number, number, number]> = [];
    const threshold = 0.5; // Å

    for (let i = 0; i < atoms.chemical_symbols.length; i++) {
      for (let j = i + 1; j < atoms.chemical_symbols.length; j++) {
        const dist = this.distance(atoms.positions[i], atoms.positions[j]);
        if (dist < threshold) {
          overlaps.push([i, j, dist]);
        }
      }
    }

    return {
      valid: overlaps.length === 0,
      overlaps,
      warnings: overlaps.map(([i, j, d]) => `Atoms ${i}-${j} overlap: ${d.toFixed(3)} Å`),
    };
  }

  /**
   * Check charge neutrality (simplified)
   */
  private checkChargeNeutrality(atoms: ASEAtoms): any {
    const charges = this.extractCharges(atoms);
    if (!charges) {
      return {
        valid: true,
        status: 'skipped',
        netCharge: null,
        warnings: ['No per-atom charges available; skipped charge neutrality check'],
      };
    }

    if (charges.length !== atoms.chemical_symbols.length) {
      return {
        valid: false,
        status: 'skipped',
        netCharge: null,
        warnings: [
          `Charge array length ${charges.length} does not match atom count ${atoms.chemical_symbols.length}`,
        ],
      };
    }

    const netCharge = charges.reduce((sum, charge) => sum + charge, 0);
    const tolerance = Number(atoms.info?.chargeTolerance ?? 1e-6);
    const valid = Math.abs(netCharge) <= tolerance;

    return {
      valid,
      status: 'checked',
      netCharge,
      tolerance,
      warnings: valid ? [] : [`Net charge is ${netCharge.toFixed(6)} e`],
    };
  }

  /**
   * Calculate center of mass
   */
  private calculateCenterOfMass(atoms: ASEAtoms): number[] {
    const n = atoms.chemical_symbols.length;
    const com = [0, 0, 0];
    const masses = this.getMasses(atoms);
    const totalMass = masses.reduce((sum, mass) => sum + mass, 0);

    if (n === 0 || totalMass === 0) {
      return com;
    }

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < 3; j++) {
        com[j] += atoms.positions[i][j] * masses[i];
      }
    }

    return com.map(v => this.cleanNumber(v / totalMass));
  }

  /**
   * Calculate moment of inertia tensor about the center of mass.
   */
  private calculateMomentOfInertia(atoms: ASEAtoms): number[][] {
    const com = this.calculateCenterOfMass(atoms);
    const masses = this.getMasses(atoms);
    const tensor = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];

    for (let i = 0; i < atoms.chemical_symbols.length; i++) {
      const mass = masses[i];
      const x = atoms.positions[i][0] - com[0];
      const y = atoms.positions[i][1] - com[1];
      const z = atoms.positions[i][2] - com[2];

      tensor[0][0] += mass * (y * y + z * z);
      tensor[1][1] += mass * (x * x + z * z);
      tensor[2][2] += mass * (x * x + y * y);
      tensor[0][1] -= mass * x * y;
      tensor[0][2] -= mass * x * z;
      tensor[1][2] -= mass * y * z;
    }

    tensor[1][0] = tensor[0][1];
    tensor[2][0] = tensor[0][2];
    tensor[2][1] = tensor[1][2];

    return tensor.map(row => row.map(value => this.cleanNumber(value)));
  }

  /**
   * Calculate bounding box
   */
  private calculateBoundingBox(atoms: ASEAtoms): any {
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];

    for (const pos of atoms.positions) {
      for (let i = 0; i < 3; i++) {
        min[i] = Math.min(min[i], pos[i]);
        max[i] = Math.max(max[i], pos[i]);
      }
    }

    return {
      min,
      max,
      size: max.map((v, i) => v - min[i]),
    };
  }

  /**
   * Calculate Euclidean distance
   */
  private distance(a: number[], b: number[]): number {
    let sum = 0;
    for (let i = 0; i < a.length; i++) {
      const diff = a[i] - b[i];
      sum += diff * diff;
    }
    return Math.sqrt(sum);
  }

  private parseXyz(content: string, format: string): ASEAtoms {
    const rawLines = content.split(/\r?\n/);
    const lines = rawLines.map(line => line.trim()).filter(line => line.length > 0);
    if (lines.length === 0) {
      throw new Error('Invalid XYZ: no atom lines found');
    }

    const declaredCount = Number.parseInt(lines[0], 10);
    const hasHeader = Number.isInteger(declaredCount) && String(declaredCount) === lines[0];
    const comment = hasHeader ? (rawLines[1] ?? '').trim() : '';
    const atomLines = hasHeader ? rawLines.slice(2, declaredCount + 2) : rawLines;
    const parsed = this.parseAtomCoordinateLines(atomLines);

    if (hasHeader && parsed.chemical_symbols.length !== declaredCount) {
      throw new Error(
        `Invalid XYZ: declared ${declaredCount} atoms but parsed ${parsed.chemical_symbols.length}`
      );
    }
    if (parsed.chemical_symbols.length === 0) {
      throw new Error('Invalid XYZ: no valid atom coordinate lines found');
    }

    return {
      chemical_symbols: parsed.chemical_symbols,
      positions: parsed.positions,
      pbc: [false, false, false],
      info: {
        format,
        comment,
        parsed: true,
        parser: 'compute-worker-native',
      },
    };
  }

  private parsePoscar(content: string, format: string): ASEAtoms {
    const lines = content
      .split(/\r?\n/)
      .map(line => line.trim())
      .filter(line => line.length > 0);
    if (lines.length < 8) {
      throw new Error('Invalid POSCAR: too few lines');
    }

    const title = lines[0];
    const scale = this.parseFiniteNumber(lines[1], 'POSCAR scale');
    const latticeScale = scale === 0 ? 1 : scale;
    const cell = [
      this.parseVector(lines[2], latticeScale, 'POSCAR cell a'),
      this.parseVector(lines[3], latticeScale, 'POSCAR cell b'),
      this.parseVector(lines[4], latticeScale, 'POSCAR cell c'),
    ];

    const line5 = lines[5].split(/\s+/);
    const line6 = lines[6].split(/\s+/);
    const hasElementLine = line5.every(token => /^[A-Za-z][A-Za-z0-9_+-]*$/.test(token));
    const elements = hasElementLine ? line5.map(element => this.normalizeElement(element)) : [];
    const counts = (hasElementLine ? line6 : line5).map((token, index) => {
      const count = Number.parseInt(token, 10);
      if (!Number.isInteger(count) || count < 0) {
        throw new Error(`Invalid POSCAR atom count at index ${index + 1}`);
      }
      return count;
    });
    const symbols = hasElementLine ? elements : counts.map((_, index) => `X${index + 1}`);

    let cursor = hasElementLine ? 7 : 6;
    let mode = lines[cursor].toLowerCase();
    if (mode.startsWith('s')) {
      cursor += 1;
      mode = lines[cursor].toLowerCase();
    }
    const direct = mode.startsWith('d');
    const cartesian = mode.startsWith('c') || mode.startsWith('k');
    if (!direct && !cartesian) {
      throw new Error(`Invalid POSCAR coordinate mode: ${lines[cursor]}`);
    }
    cursor += 1;

    const chemicalSymbols: string[] = [];
    const positions: number[][] = [];
    for (let speciesIndex = 0; speciesIndex < symbols.length; speciesIndex++) {
      for (let atomIndex = 0; atomIndex < counts[speciesIndex]; atomIndex++) {
        if (cursor >= lines.length) {
          throw new Error('Invalid POSCAR: not enough coordinate lines');
        }
        const fractionalOrCartesian = this.parseVector(lines[cursor], 1, 'POSCAR atom position');
        chemicalSymbols.push(symbols[speciesIndex]);
        positions.push(
          direct
            ? this.fractionalToCartesian(fractionalOrCartesian, cell)
            : fractionalOrCartesian.map(value => value * latticeScale)
        );
        cursor += 1;
      }
    }

    return {
      chemical_symbols: chemicalSymbols,
      positions: positions.map(position => position.map(value => this.cleanNumber(value))),
      cell,
      pbc: [true, true, true],
      info: {
        format,
        title,
        coordinateMode: direct ? 'direct' : 'cartesian',
        parsed: true,
        parser: 'compute-worker-native',
      },
    };
  }

  private formatXyz(atoms: ASEAtoms, extended: boolean): string {
    const comment =
      extended && atoms.cell
        ? `Lattice="${atoms.cell
            .flatMap(vector => vector.slice(0, 3))
            .map(value => this.formatCoordinate(value))
            .join(' ')}" Properties=species:S:1:pos:R:3 pbc="${atoms.pbc
            .map(value => (value ? 'T' : 'F'))
            .join(' ')}"`
        : String(atoms.info?.comment ?? 'OpenQC worker export');

    const lines = [String(atoms.chemical_symbols.length), comment];
    atoms.chemical_symbols.forEach((symbol, index) => {
      lines.push(
        `${this.normalizeElement(symbol)} ${atoms.positions[index]
          .slice(0, 3)
          .map(value => this.formatCoordinate(value))
          .join(' ')}`
      );
    });

    return lines.join('\n');
  }

  private formatPoscar(atoms: ASEAtoms): string {
    if (!atoms.cell) {
      throw new Error('Cannot convert to POSCAR/VASP without unit-cell vectors');
    }

    const groups = this.groupAtomIndicesBySymbol(atoms);
    const lines = [
      String(atoms.info?.title ?? 'OpenQC worker export'),
      '1.000000000000',
      ...atoms.cell.map(vector =>
        vector
          .slice(0, 3)
          .map(value => this.formatCoordinate(value))
          .join(' ')
      ),
      groups.map(group => group.symbol).join(' '),
      groups.map(group => String(group.indices.length)).join(' '),
      'Cartesian',
    ];

    groups.forEach(group => {
      group.indices.forEach(index => {
        lines.push(
          atoms.positions[index]
            .slice(0, 3)
            .map(value => this.formatCoordinate(value))
            .join(' ')
        );
      });
    });

    return lines.join('\n');
  }

  private groupAtomIndicesBySymbol(atoms: ASEAtoms): Array<{ symbol: string; indices: number[] }> {
    const groups: Array<{ symbol: string; indices: number[] }> = [];
    const bySymbol = new Map<string, number[]>();

    atoms.chemical_symbols.forEach((rawSymbol, index) => {
      const symbol = this.normalizeElement(rawSymbol);
      let indices = bySymbol.get(symbol);
      if (!indices) {
        indices = [];
        bySymbol.set(symbol, indices);
        groups.push({ symbol, indices });
      }
      indices.push(index);
    });

    return groups;
  }

  private parseAtomCoordinateLines(
    lines: string[]
  ): Pick<ASEAtoms, 'chemical_symbols' | 'positions'> {
    const chemicalSymbols: string[] = [];
    const positions: number[][] = [];
    for (const line of lines) {
      const stripped = line.trim();
      if (!stripped || stripped.startsWith('#')) {
        continue;
      }
      const parts = stripped.split(/\s+/);
      if (parts.length < 4 || !/^[A-Za-z][A-Za-z0-9_+-]*$/.test(parts[0])) {
        continue;
      }
      const position = [
        Number.parseFloat(parts[1]),
        Number.parseFloat(parts[2]),
        Number.parseFloat(parts[3]),
      ];
      if (position.some(value => !Number.isFinite(value))) {
        continue;
      }
      chemicalSymbols.push(this.normalizeElement(parts[0]));
      positions.push(position);
    }

    return { chemical_symbols: chemicalSymbols, positions };
  }

  private parseVector(line: string, scale: number, label: string): number[] {
    const parts = line.split(/\s+/);
    if (parts.length < 3) {
      throw new Error(`Invalid ${label}: expected 3 numeric values`);
    }
    return [
      this.parseFiniteNumber(parts[0], label) * scale,
      this.parseFiniteNumber(parts[1], label) * scale,
      this.parseFiniteNumber(parts[2], label) * scale,
    ];
  }

  private parseFiniteNumber(value: string, label: string): number {
    const number = Number.parseFloat(value);
    if (!Number.isFinite(number)) {
      throw new Error(`Invalid ${label}: ${value}`);
    }
    return number;
  }

  private fractionalToCartesian(fractional: number[], cell: number[][]): number[] {
    return [
      fractional[0] * cell[0][0] + fractional[1] * cell[1][0] + fractional[2] * cell[2][0],
      fractional[0] * cell[0][1] + fractional[1] * cell[1][1] + fractional[2] * cell[2][1],
      fractional[0] * cell[0][2] + fractional[1] * cell[1][2] + fractional[2] * cell[2][2],
    ];
  }

  private normalizeElement(value: string): string {
    const letters = value.replace(/[^A-Za-z]/g, '').slice(0, 2);
    if (!letters) {
      return 'X';
    }
    return letters.charAt(0).toUpperCase() + letters.slice(1).toLowerCase();
  }

  private extractCharges(atoms: ASEAtoms): number[] | undefined {
    const candidate =
      (atoms as any).charges ??
      atoms.info?.charges ??
      atoms.info?.initial_charges ??
      atoms.info?.formal_charges;
    if (!Array.isArray(candidate)) {
      return undefined;
    }

    const charges = candidate.map(value => Number(value));
    if (charges.some(value => !Number.isFinite(value))) {
      return undefined;
    }
    return charges;
  }

  private getMasses(atoms: ASEAtoms): number[] {
    if (
      Array.isArray(atoms.masses) &&
      atoms.masses.length === atoms.chemical_symbols.length &&
      atoms.masses.every(mass => Number.isFinite(mass) && mass > 0)
    ) {
      return atoms.masses;
    }

    return atoms.chemical_symbols.map(symbol => ATOMIC_MASSES[symbol] ?? 1);
  }

  private cleanNumber(value: number): number {
    return Math.abs(value) < 1e-12 ? 0 : Number(value.toFixed(12));
  }

  private cleanConvertedValue(value: any): any {
    if (typeof value === 'number' && Number.isFinite(value)) {
      return this.cleanNumber(value);
    }
    return value;
  }

  private formatCoordinate(value: number): string {
    return this.cleanNumber(value).toFixed(6);
  }

  private getBaseParameterName(key: string): string {
    const parts = key.split('.');
    return parts[parts.length - 1];
  }

  private validateAtoms(atoms: ASEAtoms | null | undefined): void {
    if (!atoms || !Array.isArray(atoms.chemical_symbols) || !Array.isArray(atoms.positions)) {
      throw new Error('Invalid atoms object');
    }
    if (atoms.chemical_symbols.length === 0) {
      throw new Error('atoms must contain at least one atom');
    }
    if (atoms.chemical_symbols.length !== atoms.positions.length) {
      throw new Error('chemical_symbols length must match positions length');
    }
    if (!Array.isArray(atoms.pbc) || atoms.pbc.length !== 3) {
      throw new Error('pbc must be a 3-element boolean array');
    }

    atoms.positions.forEach((position, index) => {
      if (
        !Array.isArray(position) ||
        position.length < 3 ||
        position.slice(0, 3).some(value => !Number.isFinite(value))
      ) {
        throw new Error(`Atom ${index}: position must contain three finite numbers`);
      }
    });

    if (atoms.cell !== undefined && atoms.cell !== null) {
      if (
        !Array.isArray(atoms.cell) ||
        atoms.cell.length !== 3 ||
        atoms.cell.some(
          vector =>
            !Array.isArray(vector) ||
            vector.length < 3 ||
            vector.slice(0, 3).some(value => !Number.isFinite(value))
        )
      ) {
        throw new Error('cell must contain three finite 3D vectors');
      }
    }
  }
}

const ATOMIC_MASSES: Record<string, number> = {
  H: 1.008,
  He: 4.002602,
  Li: 6.94,
  Be: 9.0121831,
  B: 10.81,
  C: 12.011,
  N: 14.007,
  O: 15.999,
  F: 18.998403163,
  Ne: 20.1797,
  Na: 22.98976928,
  Mg: 24.305,
  Al: 26.9815385,
  Si: 28.085,
  P: 30.973761998,
  S: 32.06,
  Cl: 35.45,
  Ar: 39.948,
  K: 39.0983,
  Ca: 40.078,
};

// Worker instance
const worker = new ComputeWorker();

// Message handler (for WebWorker environment)
if (typeof self !== 'undefined' && typeof (self as any).postMessage === 'function') {
  self.onmessage = async (event: MessageEvent) => {
    const message = event.data as WorkerMessage;
    const response = await worker.processMessage(message);
    (self as any).postMessage(response);
  };
}

// Export for testing
export { ComputeWorker };
