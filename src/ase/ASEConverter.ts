/**
 * ASE Converter - TypeScript Interface
 *
 * Provides TypeScript interface to the Python ASE backend for
 * format conversion between quantum chemistry files.
 */

import * as vscode from 'vscode';
import * as path from 'path';
import { spawn, exec } from 'child_process';
import { promisify } from 'util';
import * as fs from 'fs';

const execAsync = promisify(exec);

/**
 * Supported quantum chemistry formats
 */
export enum ASEFormat {
  VASP = 'vasp',
  CP2K = 'cp2k',
  QE = 'qe',
  Gaussian = 'gaussian',
  ORCA = 'orca',
  NWChem = 'nwchem',
  GAMESS = 'gamess',
  LAMMPS = 'lammps',
  XYZ = 'xyz',
  PDB = 'pdb',
  CIF = 'cif',
}

/**
 * ASE Atoms object representation
 */
export interface ASEAtoms {
  chemical_symbols: string[];
  positions: number[][];
  cell?: number[][] | null;
  pbc: boolean[];
  info?: Record<string, any>;
  masses?: number[];
  constraints?: string;
}

/**
 * Result of conversion operation
 */
export interface ConversionResult {
  success: boolean;
  data_type?: 'atoms' | 'text' | 'dict';
  atoms?: ASEAtoms;
  content?: string;
  data?: any;
  error?: string;
  warnings: string[];
  metadata: Record<string, any>;
}

/**
 * ASE Converter class
 *
 * Interfaces with Python ASE backend for format conversion
 */
export class ASEConverter {
  private pythonPath: string;
  private converterScript: string;
  private context: vscode.ExtensionContext;

  constructor(context: vscode.ExtensionContext) {
    this.context = context;
    this.pythonPath = this.getPythonPath();
    this.converterScript = path.join(
      context.extensionPath,
      'python',
      'openqc',
      'ase',
      'converter.py'
    );
  }

  /**
   * Get Python executable path from configuration
   */
  private getPythonPath(): string {
    const config = vscode.workspace.getConfiguration('openqc');
    return config.get<string>('pythonPath', 'python3');
  }

  /**
   * Read file to ASE Atoms
   */
  public async readToAtoms(filepath: string, formatHint?: ASEFormat): Promise<ConversionResult> {
    const args = ['read', filepath];
    if (formatHint) {
      args.push(formatHint);
    }

    return this.executeConverter(args);
  }

  /**
   * Write ASE Atoms to file
   */
  public async writeFromAtoms(
    atoms: ASEAtoms,
    outputPath: string,
    formatHint?: ASEFormat
  ): Promise<ConversionResult> {
    const args = ['write', JSON.stringify(atoms), outputPath];
    if (formatHint) {
      args.push(formatHint);
    }

    return this.executeConverter(args);
  }

  /**
   * Convert between formats
   */
  public async convertFormat(
    inputPath: string,
    outputPath: string,
    inputFormat?: ASEFormat,
    outputFormat?: ASEFormat
  ): Promise<ConversionResult> {
    const args = ['convert', inputPath, outputPath];
    if (inputFormat) {
      args.push(inputFormat);
    }
    if (outputFormat) {
      args.push(outputFormat);
    }

    return this.executeConverter(args);
  }

  /**
   * Execute Python converter script
   */
  private async executeConverter(args: string[]): Promise<ConversionResult> {
    return new Promise(resolve => {
      const process = spawn(this.pythonPath, [this.converterScript, ...args]);

      let stdout = '';
      let stderr = '';

      process.stdout.on('data', data => {
        stdout += data.toString();
      });

      process.stderr.on('data', data => {
        stderr += data.toString();
      });

      process.on('close', code => {
        if (code !== 0) {
          resolve({
            success: false,
            error: `Converter failed with code ${code}: ${stderr}`,
            warnings: [],
            metadata: {},
          });
          return;
        }

        try {
          const result = JSON.parse(stdout) as ConversionResult;
          resolve(result);
        } catch (error) {
          resolve({
            success: false,
            error: `Failed to parse converter output: ${error}`,
            warnings: [],
            metadata: { raw_output: stdout },
          });
        }
      });

      process.on('error', error => {
        resolve({
          success: false,
          error: `Failed to execute converter: ${error.message}`,
          warnings: [],
          metadata: {},
        });
      });
    });
  }

  /**
   * Check if Python backend is available
   */
  public async isAvailable(): Promise<boolean> {
    try {
      const { stdout } = await execAsync(
        `${this.pythonPath} -c "import ase; print(ase.__version__)"`
      );
      return stdout.trim().length > 0;
    } catch {
      return false;
    }
  }

  /**
   * Get supported formats
   */
  public getSupportedFormats(): Record<
    string,
    { name: string; extensions: string[]; description: string }
  > {
    return {
      vasp: {
        name: 'VASP',
        extensions: ['POSCAR', 'CONTCAR'],
        description: 'VASP structure format',
      },
      cp2k: {
        name: 'CP2K',
        extensions: ['.inp'],
        description: 'CP2K input format',
      },
      qe: {
        name: 'Quantum ESPRESSO',
        extensions: ['.in', '.pw.in'],
        description: 'Quantum ESPRESSO input format',
      },
      gaussian: {
        name: 'Gaussian',
        extensions: ['.com', '.gjf'],
        description: 'Gaussian input format',
      },
      orca: {
        name: 'ORCA',
        extensions: ['.inp'],
        description: 'ORCA input format',
      },
      nwchem: {
        name: 'NWChem',
        extensions: ['.nw', '.nwinp'],
        description: 'NWChem input format',
      },
      gamess: {
        name: 'GAMESS',
        extensions: ['.inp'],
        description: 'GAMESS input format',
      },
      lammps: {
        name: 'LAMMPS',
        extensions: ['.lmp', '.data'],
        description: 'LAMMPS data format',
      },
      xyz: {
        name: 'XYZ',
        extensions: ['.xyz'],
        description: 'Generic XYZ format',
      },
      pdb: {
        name: 'PDB',
        extensions: ['.pdb'],
        description: 'Protein Data Bank format',
      },
      cif: {
        name: 'CIF',
        extensions: ['.cif'],
        description: 'Crystallographic Information File',
      },
    };
  }
}
