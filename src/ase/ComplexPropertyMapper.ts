/**
 * ASE Complex Property Mapper
 *
 * Handles complex property mapping between quantum chemistry codes:
 * - Hubbard U parameter mapping
 * - Constraints and frozen atoms handling
 * - Excited state method mapping
 * - Pseudopotential/basis set strategies
 */

import * as vscode from 'vscode';
import * as path from 'path';
import { spawn } from 'child_process';
import { promisify } from 'util';
import { exec } from 'child_process';

const execAsync = promisify(exec);

/**
 * Hubbard U parameters
 */
export interface HubbardParameters {
  enabled: boolean;
  method: 'Dudarev' | 'Liechtenstein' | 'simplified' | 'full' | 'DUDAREV';
  uValues: Record<string, number>;
  jValues?: Record<string, number>;
}

/**
 * Constraint information
 */
export interface ConstraintInfo {
  enabled: boolean;
  constrainedAtoms: number[];
  constraintType: 'fixed' | 'line' | 'plane' | 'none';
  constraintVectors?: number[][];
}

/**
 * Excited state information
 */
export interface ExcitedStateInfo {
  enabled: boolean;
  method: 'TD-DFT' | 'GW' | 'BSE' | 'turboTDDFT' | 'TDDFPT' | 'none';
  nStates?: number;
  nOcc?: number;
  nVirt?: number;
}

/**
 * Pseudopotential information
 */
export interface PseudopotentialInfo {
  strategy: 'PAW' | 'NORMCONSERVING' | 'GTH';
  functional: string;
  elementMapping: Record<string, string>;
  basisSetMapping?: Record<string, string>;
}

/**
 * Supported quantum chemistry codes
 */
export type QuantumCode = 'vasp' | 'qe' | 'cp2k' | 'gaussian' | 'orca';

/**
 * Complex Property Mapper class
 *
 * Maps complex properties between quantum chemistry codes
 */
export class ComplexPropertyMapper {
  private pythonPath: string;
  private mapperScript: string;
  private context: vscode.ExtensionContext;

  constructor(context: vscode.ExtensionContext) {
    this.context = context;
    this.pythonPath = this.getPythonPath();
    this.mapperScript = path.join(
      context.extensionPath,
      'python',
      'openqc',
      'ase',
      'complex_properties.py'
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
   * Convert Hubbard U parameters between codes
   */
  public async convertHubbard(
    hubbard: HubbardParameters,
    sourceCode: QuantumCode,
    targetCode: QuantumCode
  ): Promise<HubbardParameters> {
    const args = [
      'hubbard',
      sourceCode,
      targetCode,
      JSON.stringify({
        enabled: hubbard.enabled,
        method: hubbard.method,
        u_values: hubbard.uValues,
        j_values: hubbard.jValues,
      }),
    ];

    const result = await this.executeMapper(args);
    if (!result.success) {
      throw new Error(result.error || 'Failed to convert Hubbard parameters');
    }

    return {
      enabled: result.data.enabled,
      method: result.data.method,
      uValues: result.data.u_values,
      jValues: result.data.j_values,
    };
  }

  /**
   * Convert constraints between codes
   */
  public async convertConstraints(
    constraints: ConstraintInfo,
    sourceCode: QuantumCode,
    targetCode: QuantumCode
  ): Promise<ConstraintInfo> {
    const args = [
      'constraints',
      sourceCode,
      targetCode,
      JSON.stringify({
        enabled: constraints.enabled,
        constrained_atoms: constraints.constrainedAtoms,
        constraint_type: constraints.constraintType,
        constraint_vectors: constraints.constraintVectors,
      }),
    ];

    const result = await this.executeMapper(args);
    if (!result.success) {
      throw new Error(result.error || 'Failed to convert constraints');
    }

    return {
      enabled: result.data.enabled,
      constrainedAtoms: result.data.constrained_atoms,
      constraintType: result.data.constraint_type,
      constraintVectors: result.data.constraint_vectors,
    };
  }

  /**
   * Convert excited state method between codes
   */
  public async convertExcitedState(
    excited: ExcitedStateInfo,
    sourceCode: QuantumCode,
    targetCode: QuantumCode
  ): Promise<ExcitedStateInfo> {
    const args = [
      'excited',
      sourceCode,
      targetCode,
      JSON.stringify({
        enabled: excited.enabled,
        method: excited.method,
        n_states: excited.nStates,
        n_occ: excited.nOcc,
        n_virt: excited.nVirt,
      }),
    ];

    const result = await this.executeMapper(args);
    if (!result.success) {
      throw new Error(result.error || 'Failed to convert excited state method');
    }

    return {
      enabled: result.data.enabled,
      method: result.data.method,
      nStates: result.data.n_states,
      nOcc: result.data.n_occ,
      nVirt: result.data.n_virt,
    };
  }

  /**
   * Convert pseudopotential information between codes
   */
  public async convertPseudopotential(
    ppInfo: PseudopotentialInfo,
    sourceCode: QuantumCode,
    targetCode: QuantumCode
  ): Promise<PseudopotentialInfo> {
    const args = [
      'pp',
      sourceCode,
      targetCode,
      JSON.stringify({
        strategy: ppInfo.strategy,
        functional: ppInfo.functional,
        element_mapping: ppInfo.elementMapping,
        basis_set_mapping: ppInfo.basisSetMapping,
      }),
    ];

    const result = await this.executeMapper(args);
    if (!result.success) {
      throw new Error(result.error || 'Failed to convert pseudopotential');
    }

    return {
      strategy: result.data.strategy,
      functional: result.data.functional,
      elementMapping: result.data.element_mapping,
      basisSetMapping: result.data.basis_set_mapping,
    };
  }

  /**
   * Generate VASP Hubbard section for INCAR
   */
  public async generateVASPHubbardSection(hubbard: HubbardParameters): Promise<string> {
    const args = [
      'generate_vasp_hubbard',
      JSON.stringify({
        enabled: hubbard.enabled,
        method: hubbard.method,
        u_values: hubbard.uValues,
        j_values: hubbard.jValues,
      }),
    ];

    const result = await this.executeMapper(args);
    return result.data || '';
  }

  /**
   * Generate QE Hubbard section for input
   */
  public async generateQEHubbardSection(hubbard: HubbardParameters): Promise<string> {
    const args = [
      'generate_qe_hubbard',
      JSON.stringify({
        enabled: hubbard.enabled,
        method: hubbard.method,
        u_values: hubbard.uValues,
        j_values: hubbard.jValues,
      }),
    ];

    const result = await this.executeMapper(args);
    return result.data || '';
  }

  /**
   * Generate CP2K Hubbard section for input
   */
  public async generateCP2KHubbardSection(hubbard: HubbardParameters): Promise<string> {
    const args = [
      'generate_cp2k_hubbard',
      JSON.stringify({
        enabled: hubbard.enabled,
        method: hubbard.method,
        u_values: hubbard.uValues,
        j_values: hubbard.jValues,
      }),
    ];

    const result = await this.executeMapper(args);
    return result.data || '';
  }

  /**
   * Check if mapper is available
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
   * Execute Python mapper script
   */
  private async executeMapper(args: string[]): Promise<{ success: boolean; data?: any; error?: string }> {
    return new Promise((resolve) => {
      const process = spawn(this.pythonPath, [this.mapperScript, ...args]);

      let stdout = '';
      let stderr = '';

      process.stdout.on('data', (data) => {
        stdout += data.toString();
      });

      process.stderr.on('data', (data) => {
        stderr += data.toString();
      });

      process.on('close', (code) => {
        if (code !== 0) {
          resolve({
            success: false,
            error: `Mapper failed with code ${code}: ${stderr}`,
          });
          return;
        }

        try {
          const data = JSON.parse(stdout);
          resolve({ success: true, data });
        } catch {
          // For generate commands, return raw output
          resolve({ success: true, data: stdout });
        }
      });

      process.on('error', (error) => {
        resolve({
          success: false,
          error: `Failed to execute mapper: ${error.message}`,
        });
      });
    });
  }
}

/**
 * Export singleton instance
 */
export let complexPropertyMapper: ComplexPropertyMapper;

/**
 * Initialize complex property mapper
 */
export function initializeComplexPropertyMapper(
  context: vscode.ExtensionContext
): ComplexPropertyMapper {
  complexPropertyMapper = new ComplexPropertyMapper(context);
  return complexPropertyMapper;
}
