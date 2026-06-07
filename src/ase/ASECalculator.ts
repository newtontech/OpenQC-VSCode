/**
 * ASE Calculator - TypeScript Interface
 *
 * Provides TypeScript interface to ASE calculators (VASP, CP2K, QE)
 * for generating input files, executing calculations, and reading results.
 */

import * as vscode from 'vscode';
import * as path from 'path';
import { spawn, exec } from 'child_process';
import { promisify } from 'util';
import * as fs from 'fs';
import { ASEFormat, ASEAtoms, ConversionResult } from './ASEConverter';

const execAsync = promisify(exec);

/**
 * Supported calculator types
 */
export enum CalculatorType {
  VASP = 'vasp',
  CP2K = 'cp2k',
  QE = 'qe',
}

/**
 * Calculator configuration interface
 */
export interface CalculatorConfig {
  type: CalculatorType;
  command?: string;
  directory?: string;
  environment?: Record<string, string>;
  parameters: Record<string, any>;
  pseudopotentials?: Record<string, string>;
  basisSets?: Record<string, string>;
}

/**
 * Calculation result interface
 */
export interface CalculationResult {
  success: boolean;
  energy?: number;
  forces?: number[][];
  stress?: number[][];
  atoms?: ASEAtoms;
  error?: string;
  warnings: string[];
  metadata: Record<string, any>;
  outputFiles: string[];
}

/**
 * Input generation result
 */
export interface InputGenerationResult {
  success: boolean;
  files: Record<string, string>;
  error?: string;
  warnings: string[];
}

/**
 * ASE Calculator class
 *
 * Interfaces with Python ASE backend for calculator operations
 */
export class ASECalculator {
  private pythonPath: string;
  private calculatorScript: string;
  private context: vscode.ExtensionContext;
  private config: CalculatorConfig;

  constructor(context: vscode.ExtensionContext, config: CalculatorConfig) {
    this.context = context;
    this.config = config;
    this.pythonPath = this.getPythonPath();
    this.calculatorScript = path.join(
      context.extensionPath,
      'python',
      'openqc',
      'ase',
      'calculator.py'
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
   * Generate input files for a calculation
   */
  public async generateInput(atoms: ASEAtoms, outputDir: string): Promise<InputGenerationResult> {
    const args = [
      'generate',
      this.config.type,
      JSON.stringify(atoms),
      outputDir,
      JSON.stringify(this.config.parameters),
    ];

    if (this.config.pseudopotentials) {
      args.push('--pseudos', JSON.stringify(this.config.pseudopotentials));
    }

    if (this.config.basisSets) {
      args.push('--basis', JSON.stringify(this.config.basisSets));
    }

    return this.executeCalculator<InputGenerationResult>(args);
  }

  /**
   * Run a calculation
   */
  public async runCalculation(inputDir: string, outputDir?: string): Promise<CalculationResult> {
    const args = ['run', this.config.type, inputDir];

    if (outputDir) {
      args.push('--output-dir', outputDir);
    }

    if (this.config.command) {
      args.push('--command', this.config.command);
    }

    if (this.config.environment) {
      args.push('--env', JSON.stringify(this.config.environment));
    }

    return this.executeCalculator<CalculationResult>(args);
  }

  /**
   * Read calculation results
   */
  public async readResults(outputDir: string): Promise<CalculationResult> {
    const args = ['read', this.config.type, outputDir];
    return this.executeCalculator<CalculationResult>(args);
  }

  /**
   * Generate input and run calculation in one step
   */
  public async calculate(atoms: ASEAtoms, workingDir: string): Promise<CalculationResult> {
    // Generate input files
    const inputResult = await this.generateInput(atoms, workingDir);
    if (!inputResult.success) {
      return {
        success: false,
        error: `Input generation failed: ${inputResult.error}`,
        warnings: inputResult.warnings,
        metadata: {},
        outputFiles: [],
      };
    }

    // Run calculation
    const runResult = await this.runCalculation(workingDir);
    if (!runResult.success) {
      return {
        success: false,
        error: `Calculation failed: ${runResult.error}`,
        warnings: [...inputResult.warnings, ...runResult.warnings],
        metadata: runResult.metadata,
        outputFiles: runResult.outputFiles,
      };
    }

    return runResult;
  }

  /**
   * Execute Python calculator script
   */
  private async executeCalculator<T>(args: string[]): Promise<T> {
    return new Promise(resolve => {
      const process = spawn(this.pythonPath, [this.calculatorScript, ...args]);

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
            error: `Calculator failed with code ${code}: ${stderr}`,
            warnings: [],
            outputFiles: [],
            metadata: {},
          } as unknown as T);
          return;
        }

        try {
          const result = JSON.parse(stdout) as T;
          resolve(result);
        } catch (error) {
          resolve({
            success: false,
            error: `Failed to parse calculator output: ${error}`,
            warnings: [],
            outputFiles: [],
            metadata: { raw_output: stdout },
          } as unknown as T);
        }
      });

      process.on('error', error => {
        resolve({
          success: false,
          error: `Failed to execute calculator: ${error.message}`,
          warnings: [],
          outputFiles: [],
          metadata: {},
        } as unknown as T);
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
   * Check if calculator executable is available
   */
  public async isCalculatorAvailable(): Promise<boolean> {
    if (!this.config.command) {
      return false;
    }

    try {
      const { stdout } = await execAsync(`which ${this.config.command.split(' ')[0]}`);
      return stdout.trim().length > 0;
    } catch {
      return false;
    }
  }

  /**
   * Get default parameters for calculator type
   */
  public getDefaultParameters(): Record<string, any> {
    const defaults: Record<CalculatorType, Record<string, any>> = {
      [CalculatorType.VASP]: {
        encut: 400,
        kpts: [4, 4, 4],
        prec: 'Normal',
        ismear: 0,
        sigma: 0.1,
        ibrion: 2,
        nsw: 100,
        ediff: 1e-5,
        ediffg: -0.01,
        isif: 2,
      },
      [CalculatorType.CP2K]: {
        basis_set: 'DZVP-MOLOPT-SR-GTH',
        potential: 'GTH-PBE',
        xc_functional: 'PBE',
        cutoff: 400,
        rel_cutoff: 50,
        ngrids: 4,
      },
      [CalculatorType.QE]: {
        calculation: 'scf',
        ecutwfc: 40,
        ecutrho: 320,
        conv_thr: 1e-8,
        mixing_beta: 0.7,
        kpts: [4, 4, 4],
      },
    };

    return defaults[this.config.type] || {};
  }

  /**
   * Get calculator documentation URL
   */
  public getDocumentationUrl(): string {
    const docs: Record<CalculatorType, string> = {
      [CalculatorType.VASP]: 'https://www.vasp.at/wiki/index.php/The_VASP_Manual',
      [CalculatorType.CP2K]: 'https://manual.cp2k.org/',
      [CalculatorType.QE]: 'https://www.quantum-espresso.org/documentation/',
    };

    return docs[this.config.type] || '';
  }

  /**
   * Update calculator configuration
   */
  public updateConfig(config: Partial<CalculatorConfig>): void {
    this.config = { ...this.config, ...config };
  }

  /**
   * Get current configuration
   */
  public getConfig(): CalculatorConfig {
    return { ...this.config };
  }
}

/**
 * Calculator factory for creating calculator instances
 */
export class CalculatorFactory {
  private context: vscode.ExtensionContext;

  constructor(context: vscode.ExtensionContext) {
    this.context = context;
  }

  /**
   * Create a calculator instance
   */
  createCalculator(config: CalculatorConfig): ASECalculator {
    return new ASECalculator(this.context, config);
  }

  /**
   * Create a VASP calculator
   */
  createVASPCalculator(parameters?: Record<string, any>, command?: string): ASECalculator {
    return this.createCalculator({
      type: CalculatorType.VASP,
      command: command || 'vasp',
      parameters: parameters || {},
    });
  }

  /**
   * Create a CP2K calculator
   */
  createCP2KCalculator(parameters?: Record<string, any>, command?: string): ASECalculator {
    return this.createCalculator({
      type: CalculatorType.CP2K,
      command: command || 'cp2k.popt',
      parameters: parameters || {},
    });
  }

  /**
   * Create a Quantum ESPRESSO calculator
   */
  createQECalculator(parameters?: Record<string, any>, command?: string): ASECalculator {
    return this.createCalculator({
      type: CalculatorType.QE,
      command: command || 'pw.x',
      parameters: parameters || {},
    });
  }
}

/**
 * Get VSCode configuration for calculators
 */
export function getCalculatorConfiguration(): Record<string, any> {
  const config = vscode.workspace.getConfiguration('openqc.calculators');
  return {
    vasp: {
      command: config.get<string>('vasp.command', 'vasp'),
      defaultParams: config.get<Record<string, any>>('vasp.defaultParams', {}),
    },
    cp2k: {
      command: config.get<string>('cp2k.command', 'cp2k.popt'),
      defaultParams: config.get<Record<string, any>>('cp2k.defaultParams', {}),
    },
    qe: {
      command: config.get<string>('qe.command', 'pw.x'),
      defaultParams: config.get<Record<string, any>>('qe.defaultParams', {}),
    },
  };
}
