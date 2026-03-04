/**
 * Compute Worker - WebWorker for Heavy Computations
 *
 * Offloads computationally intensive operations to background threads
 * to prevent UI blocking in VSCode.
 */

import { ASEAtoms } from '../ase/ASEConverter';

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
    // Simulate heavy parsing (actual implementation would call Python backend)
    const { content, format, options } = payload;
    
    // Basic validation
    if (!content || content.trim().length === 0) {
      throw new Error('Empty content provided');
    }

    // For now, return a mock structure
    // In production, this would call the Python converter
    const lines = content.split('\n');
    const atomCount = Math.min(lines.length, 1000);

    const atoms: ASEAtoms = {
      chemical_symbols: Array(atomCount).fill('C'),
      positions: Array(atomCount)
        .fill(0)
        .map((_, i) => [i * 1.5, i * 1.5, i * 1.5]),
      pbc: [true, true, true],
      cell: [
        [10, 0, 0],
        [0, 10, 0],
        [0, 0, 10],
      ],
      info: { format, parsed: true },
    };

    return atoms;
  }

  /**
   * Convert atoms to target format
   */
  private async convertFormat(payload: ConvertFormatPayload): Promise<string> {
    const { atoms, targetFormat, options } = payload;

    // Validate input
    if (!atoms || !atoms.chemical_symbols || atoms.chemical_symbols.length === 0) {
      throw new Error('Invalid atoms object');
    }

    // Simulate conversion (would call Python backend)
    const lines: string[] = [];
    
    // Header
    lines.push(`# Converted to ${targetFormat}`);
    lines.push(`# Atoms: ${atoms.chemical_symbols.length}`);
    lines.push('');

    // Cell (if periodic)
    if (atoms.cell && atoms.pbc.some((p) => p)) {
      lines.push('# Unit cell:');
      atoms.cell.forEach((row) => {
        lines.push(row.map((v) => v.toFixed(6)).join(' '));
      });
      lines.push('');
    }

    // Atoms
    lines.push('# Atomic positions:');
    atoms.chemical_symbols.forEach((symbol, i) => {
      const pos = atoms.positions[i];
      lines.push(`${symbol} ${pos.map((v) => v.toFixed(6)).join(' ')}`);
    });

    return lines.join('\n');
  }

  /**
   * Validate structure
   */
  private async validateStructure(payload: ValidateStructurePayload): Promise<Record<string, any>> {
    const { atoms, checks } = payload;
    const results: Record<string, any> = {};

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
  private async calculateProperties(payload: CalculatePropertiesPayload): Promise<Record<string, any>> {
    const { atoms, properties } = payload;
    const results: Record<string, any> = {};

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

    // Simulate parameter mapping
    const migrated: Record<string, any> = {};

    // Example mappings
    if (sourceFormat === 'vasp' && targetFormat === 'cp2k') {
      if (parameters.ENCUT) {
        migrated.CUTOFF = Math.round(parameters.ENCUT * 1.1); // CP2K needs higher cutoff
      }
      if (parameters.EDIFF) {
        migrated.EPSCF = parameters.EDIFF;
      }
    } else if (sourceFormat === 'qe' && targetFormat === 'vasp') {
      if (parameters.ecutwfc) {
        migrated.ENCUT = parameters.ecutwfc;
      }
    }

    return {
      migrated,
      warnings: [],
      notes: [`Migrated from ${sourceFormat} to ${targetFormat}`],
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
    if (!atoms.cell || !atoms.pbc.some((p) => p)) {
      return { valid: true, warnings: ['Non-periodic system'] };
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
    // Simplified check - would need actual charges in production
    return {
      valid: true,
      warnings: ['Charge neutrality check not implemented'],
    };
  }

  /**
   * Calculate center of mass
   */
  private calculateCenterOfMass(atoms: ASEAtoms): number[] {
    const n = atoms.chemical_symbols.length;
    const com = [0, 0, 0];

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < 3; j++) {
        com[j] += atoms.positions[i][j];
      }
    }

    return com.map((v) => v / n);
  }

  /**
   * Calculate moment of inertia tensor (simplified)
   */
  private calculateMomentOfInertia(atoms: ASEAtoms): number[][] {
    // Simplified - return identity matrix
    return [
      [1, 0, 0],
      [0, 1, 0],
      [0, 0, 1],
    ];
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
}

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
