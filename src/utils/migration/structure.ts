/**
 * Structure Migration Tool
 *
 * Provides structure migration between quantum chemistry codes with parameter mapping
 * and validation capabilities.
 */

import * as vscode from 'vscode';
import * as path from 'path';
import { ASEConverter, ASEFormat, ASEAtoms } from '../../ase/ASEConverter';
import {
  getParameterMappings,
  convertParameterValue,
  isMigrationSupported,
  ParameterMapping,
} from './params';

/**
 * Migration result
 */
export interface MigrationResult {
  success: boolean;
  sourcePath?: string;
  targetPath?: string;
  sourceFormat?: string;
  targetFormat?: string;
  warnings: string[];
  errors: string[];
  metadata: {
    natoms: number;
    formula?: string;
    conversions: string[];
  };
}

/**
 * Migration options
 */
export interface MigrationOptions {
  /** Source format (auto-detected if not specified) */
  sourceFormat?: ASEFormat;
  /** Target format */
  targetFormat: ASEFormat;
  /** Output path (auto-generated if not specified) */
  outputPath?: string;
  /** Preserve constraints */
  preserveConstraints?: boolean;
  /** Validate output */
  validate?: boolean;
  /** Overwrite existing file */
  overwrite?: boolean;
}

/**
 * Structure Migration class
 */
export class StructureMigration {
  private converter: ASEConverter;
  private context: vscode.ExtensionContext;

  constructor(context: vscode.ExtensionContext) {
    this.context = context;
    this.converter = new ASEConverter(context);
  }

  /**
   * Migrate structure from source to target format
   */
  public async migrate(
    sourcePath: string,
    options: MigrationOptions
  ): Promise<MigrationResult> {
    const result: MigrationResult = {
      success: false,
      warnings: [],
      errors: [],
      metadata: {
        natoms: 0,
        conversions: [],
      },
    };

    try {
      // Validate migration support
      const sourceFormat = options.sourceFormat || this.detectFormat(sourcePath);
      const targetFormat = options.targetFormat;

      if (!isMigrationSupported(sourceFormat, targetFormat)) {
        result.errors.push(
          `Migration from ${sourceFormat} to ${targetFormat} is not supported`
        );
        return result;
      }

      // Generate output path if not provided
      const outputPath =
        options.outputPath ||
        this.generateOutputPath(sourcePath, targetFormat);

      result.sourcePath = sourcePath;
      result.targetPath = outputPath;
      result.sourceFormat = sourceFormat;
      result.targetFormat = targetFormat;

      // Read source structure
      const readResult = await this.converter.readToAtoms(
        sourcePath,
        sourceFormat as ASEFormat
      );

      if (!readResult.success) {
        result.errors.push(
          `Failed to read source file: ${readResult.error}`
        );
        return result;
      }

      if (!readResult.atoms) {
        result.errors.push('No atoms data found in source file');
        return result;
      }

      // Extract metadata
      result.metadata.natoms = readResult.atoms.chemical_symbols.length;
      result.metadata.formula = readResult.atoms.info?.formula;

      // Validate structure preservation
      if (options.validate !== false) {
        const validationWarnings = this.validateStructure(
          readResult.atoms,
          sourceFormat
        );
        result.warnings.push(...validationWarnings);
      }

      // Handle constraints
      let atoms = readResult.atoms;
      if (!options.preserveConstraints) {
        delete atoms.constraints;
        result.warnings.push(
          'Constraints not preserved in migration (may need manual adjustment)'
        );
      }

      // Write target structure
      const writeResult = await this.converter.writeFromAtoms(
        atoms,
        outputPath,
        targetFormat
      );

      if (!writeResult.success) {
        result.errors.push(
          `Failed to write target file: ${writeResult.error}`
        );
        return result;
      }

      // Record conversion
      result.metadata.conversions.push(
        `${sourceFormat} -> ${targetFormat}`
      );

      // Validate output if requested
      if (options.validate !== false) {
        const validation = await this.validateMigration(
          sourcePath,
          outputPath,
          sourceFormat,
          targetFormat
        );

        result.warnings.push(...validation.warnings);
        
        if (!validation.success) {
          result.errors.push(...validation.errors);
          return result;
        }
      }

      result.success = true;
      return result;
    } catch (error: any) {
      result.errors.push(`Migration failed: ${error.message}`);
      return result;
    }
  }

  /**
   * Detect format from file path
   */
  private detectFormat(filepath: string): string {
    const ext = path.extname(filepath).toLowerCase();
    const filename = path.basename(filepath).toUpperCase();

    // Special cases
    if (filename === 'POSCAR' || filename === 'CONTCAR') {
      return 'vasp';
    }

    // Extension mapping
    const formatMap: Record<string, string> = {
      '.com': 'gaussian',
      '.gjf': 'gaussian',
      '.inp': 'cp2k', // Default for .inp
      '.nw': 'nwchem',
      '.nwinp': 'nwchem',
      '.lmp': 'lammps',
      '.data': 'lammps',
      '.xyz': 'xyz',
      '.pdb': 'pdb',
      '.cif': 'cif',
    };

    return formatMap[ext] || path.basename(filepath, ext).toLowerCase();
  }

  /**
   * Generate output path
   */
  private generateOutputPath(
    sourcePath: string,
    targetFormat: string
  ): string {
    const dir = path.dirname(sourcePath);
    const baseName = path.basename(sourcePath, path.extname(sourcePath));
    const extMap: Record<string, string> = {
      vasp: 'POSCAR',
      cp2k: '.inp',
      qe: '.in',
      gaussian: '.gjf',
      orca: '.inp',
      nwchem: '.nw',
      gamess: '.inp',
      lammps: '.data',
      xyz: '.xyz',
      pdb: '.pdb',
      cif: '.cif',
    };

    const ext = extMap[targetFormat] || '.dat';
    return path.join(dir, `${baseName}_${targetFormat}${ext}`);
  }

  /**
   * Validate structure before migration
   */
  private validateStructure(atoms: ASEAtoms, format: string): string[] {
    const warnings: string[] = [];

    // Check for valid positions
    if (!atoms.positions || atoms.positions.length === 0) {
      warnings.push('No atomic positions found');
    }

    // Check cell for periodic systems
    if (atoms.pbc && atoms.pbc.some(p => p)) {
      if (!atoms.cell || atoms.cell.length !== 3) {
        warnings.push(
          'Periodic boundary conditions enabled but no valid cell found'
        );
      }
    }

    // Check for NaN or infinite values
    atoms.positions.forEach((pos, idx) => {
      pos.forEach((val, coordIdx) => {
        if (
          isNaN(val) ||
          !isFinite(val) ||
          Math.abs(val) > 1e6
        ) {
          warnings.push(
            `Invalid coordinate detected for atom ${idx + 1}: ${val}`
          );
        }
      });
    });

    return warnings;
  }

  /**
   * Validate migration by comparing structures
   */
  private async validateMigration(
    sourcePath: string,
    targetPath: string,
    sourceFormat: string,
    targetFormat: string
  ): Promise<{
    success: boolean;
    warnings: string[];
    errors: string[];
  }> {
    const result = {
      success: true,
      warnings: [],
      errors: [],
    };

    try {
      // Read both structures
      const sourceResult = await this.converter.readToAtoms(
        sourcePath,
        sourceFormat as ASEFormat
      );
      const targetResult = await this.converter.readToAtoms(
        targetPath,
        targetFormat as ASEFormat
      );

      if (!sourceResult.success || !sourceResult.atoms) {
        result.errors.push('Cannot read source file for validation');
        return result;
      }

      if (!targetResult.success || !targetResult.atoms) {
        result.errors.push('Cannot read target file for validation');
        return result;
      }

      // Compare atom counts
      const sourceNatoms = sourceResult.atoms.chemical_symbols.length;
      const targetNatoms = targetResult.atoms.chemical_symbols.length;

      if (sourceNatoms !== targetNatoms) {
        result.errors.push(
          `Atom count mismatch: source=${sourceNatoms}, target=${targetNatoms}`
        );
      }

      // Compare element composition
      const sourceElements = new Set(sourceResult.atoms.chemical_symbols);
      const targetElements = new Set(targetResult.atoms.chemical_symbols);

      if (
        sourceElements.size !== targetElements.size ||
        [...sourceElements].some(el => !targetElements.has(el))
      ) {
        result.warnings.push('Element composition may have changed');
      }

      // Compare cell if periodic
      if (sourceResult.atoms.pbc && sourceResult.atoms.pbc.some(p => p)) {
        const sourceCell = sourceResult.atoms.cell;
        const targetCell = targetResult.atoms.cell;

        if (sourceCell && targetCell) {
          // Check if cells are similar (within tolerance)
          const tolerance = 1e-4;
          for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
              if (
                Math.abs(sourceCell[i][j] - targetCell[i][j]) >
                tolerance
              ) {
                result.warnings.push(
                  `Cell vectors differ significantly at [${i}][${j}]`
                );
              }
            }
          }
        }
      }

      result.success = result.errors.length === 0;
      return result;
    } catch (error: any) {
      result.errors.push(`Validation failed: ${error.message}`);
      result.success = false;
      return result;
    }
  }

  /**
   * Get parameter mappings for a migration
   */
  public getParameterMappings(
    sourceFormat: string,
    targetFormat: string
  ): ParameterMapping[] {
    return getParameterMappings(sourceFormat, targetFormat);
  }

  /**
   * Check if migration is supported
   */
  public isSupported(sourceFormat: string, targetFormat: string): boolean {
    return isMigrationSupported(sourceFormat, targetFormat);
  }
}

/**
 * Quick migrate structure using Quick Pick
 */
export async function quickMigrateStructure(
  context: vscode.ExtensionContext
): Promise<void> {
  // Get active editor
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const sourcePath = editor.document.uri.fsPath;
  const migration = new StructureMigration(context);

  // Select target format
  const targetFormat = await vscode.window.showQuickPick(
    [
      'vasp',
      'cp2k',
      'qe',
      'gaussian',
      'orca',
      'nwchem',
      'gamess',
      'lammps',
      'xyz',
      'pdb',
      'cif',
    ].map(format => ({
      label: format.toUpperCase(),
      description: `Convert to ${format.toUpperCase()} format`,
      value: format,
    })),
    {
      placeHolder: 'Select target format',
    }
  );

  if (!targetFormat) {
    return;
  }

  // Show progress
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Migrating to ${targetFormat.label}...`,
      cancellable: false,
    },
    async () => {
      const result = await migration.migrate(sourcePath, {
        targetFormat: targetFormat.value as ASEFormat,
        validate: true,
        preserveConstraints: true,
      });

      if (result.success) {
        vscode.window.showInformationMessage(
          `Migration successful! Output: ${result.targetPath}`,
          'Open File'
        ).then(selection => {
          if (selection === 'Open File' && result.targetPath) {
            vscode.workspace.openTextDocument(result.targetPath);
          }
        });
      } else {
        vscode.window.showErrorMessage(
          `Migration failed: ${result.errors.join(', ')}`
        );
      }
    }
  );
}
