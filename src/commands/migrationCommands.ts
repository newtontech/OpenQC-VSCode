import * as path from 'path';
/**
 * Migration Commands - VSCode Extension Commands
 *
 * Provides VSCode commands for structure migration, k-point migration,
 * parameter migration, and MD workflow migration.
 */

import * as vscode from 'vscode';
import { StructureMigration, quickMigrateStructure } from '../utils/migration/structure';
import { KpointMigration, quickMigrateKpoints } from '../utils/migration/kpoints';
import { ASEFormat } from '../ase/ASEConverter';

/**
 * Register migration commands
 */
export function registerMigrationCommands(
  context: vscode.ExtensionContext
): void {
  // Structure migration command
  const migrateStructureCommand = vscode.commands.registerCommand(
    'OpenQC.migrateStructure',
    async () => {
      await quickMigrateStructure(context);
    }
  );
  context.subscriptions.push(migrateStructureCommand);

  // Advanced structure migration with options
  const migrateStructureAdvancedCommand = vscode.commands.registerCommand(
    'OpenQC.migrateStructureAdvanced',
    async () => {
      await migrateStructureWithOptions(context);
    }
  );
  context.subscriptions.push(migrateStructureAdvancedCommand);

  // k-Point grid migration
  const migrateKpointsCommand = vscode.commands.registerCommand(
    'OpenQC.migrateKpoints',
    async () => {
      await migrateKpoints(context);
    }
  );
  context.subscriptions.push(migrateKpointsCommand);

  // Electronic parameter migration
  const migrateParametersCommand = vscode.commands.registerCommand(
    'OpenQC.migrateParameters',
    async () => {
      await migrateParameters(context);
    }
  );
  context.subscriptions.push(migrateParametersCommand);

  // MD parameters migration
  const migrateMDParametersCommand = vscode.commands.registerCommand(
    'OpenQC.migrateMDParameters',
    async () => {
      await migrateMDParameters(context);
    }
  );
  context.subscriptions.push(migrateMDParametersCommand);

  // Show migration validation report
  const showValidationCommand = vscode.commands.registerCommand(
    'OpenQC.showMigrationValidation',
    async () => {
      await showMigrationValidation(context);
    }
  );
  context.subscriptions.push(showValidationCommand);
}

/**
 * Migrate structure with advanced options
 */
async function migrateStructureWithOptions(
  context: vscode.ExtensionContext
): Promise<void> {
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

  // Show options
  const preserveConstraints = await vscode.window.showQuickPick(
    [
      { label: 'Yes, preserve constraints', value: true },
      { label: 'No, remove constraints', value: false },
    ],
    {
      placeHolder: 'Preserve structural constraints?',
    }
  );

  if (preserveConstraints === undefined) {
    return;
  }

  const validate = await vscode.window.showQuickPick(
    [
      { label: 'Yes, validate output', value: true },
      { label: 'No, skip validation', value: false },
    ],
    {
      placeHolder: 'Validate migration output?',
    }
  );

  if (validate === undefined) {
    return;
  }

  // Select output location
  const outputLocation = await vscode.window.showQuickPick(
    [
      { label: 'Same directory as source', value: 'same' },
      { label: 'Browse for location', value: 'browse' },
    ],
    {
      placeHolder: 'Select output location',
    }
  );

  let outputPath: string | undefined;
  if (outputLocation?.value === 'browse') {
    const uri = await vscode.window.showSaveDialog({
      defaultUri: vscode.Uri.file(sourcePath),
      filters: {
        [targetFormat.label]: targetFormat.value === 'vasp' 
          ? ['POSCAR', 'CONTCAR'] 
          : [`*${path.extname(sourcePath)}`],
      },
    });
    if (uri) {
      outputPath = uri.fsPath;
    }
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
        outputPath,
        validate: validate.value,
        preserveConstraints: preserveConstraints.value,
      });

      if (result.success) {
        vscode.window.showInformationMessage(
          `Migration successful!\nOutput: ${result.targetPath}\nAtoms: ${result.metadata.natoms}${result.metadata.formula ? `\nFormula: ${result.metadata.formula}` : ''}`,
          'Open File',
          'Copy Path'
        ).then(selection => {
          if (selection === 'Open File' && result.targetPath) {
            vscode.workspace.openTextDocument(result.targetPath);
          } else if (selection === 'Copy Path' && result.targetPath) {
            vscode.env.clipboard.writeText(result.targetPath);
          }
        });
      } else {
        vscode.window.showErrorMessage(
          `Migration failed:\n${result.errors.join('\n')}`
        );
      }
    }
  );
}

/**
 * Migrate k-point grid
 */
async function migrateKpoints(context: vscode.ExtensionContext): Promise<void> {
  await quickMigrateKpoints(context);
}
/**
 * Migrate electronic structure parameters
 */
async function migrateParameters(context: vscode.ExtensionContext): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const sourcePath = editor.document.uri.fsPath;
  const migration = new StructureMigration(context);

  // Detect format from file
  const ext = sourcePath.split('.').pop()?.toLowerCase() || '';
  const formatMap: Record<string, string> = {
    'incar': 'vasp',
    'inp': 'cp2k',
    'in': 'qe',
    'com': 'gaussian',
    'gjf': 'gaussian',
    'nw': 'nwchem',
  };
  const sourceFormat = formatMap[ext] || ext;

  // Select target format
  const targetFormat = await vscode.window.showQuickPick(
    ['vasp', 'cp2k', 'qe', 'gaussian'].map(code => ({
      label: code.toUpperCase(),
      value: code,
    })),
    {
      placeHolder: 'Select target format',
    }
  );

  if (!targetFormat) {
    return;
  }

  // Get parameter mappings
  const mappings = migration.getParameterMappings(
    sourceFormat,
    targetFormat.value
  );

  if (mappings.length === 0) {
    vscode.window.showInformationMessage(
      `No parameter mappings found for ${sourceFormat} -> ${targetFormat.label}`
    );
    return;
  }

  // Show mappings
  const mappingsText = mappings
    .map(
      m =>
        `${m.sourceParam} -> ${m.targetParam}\n  ${m.description}`
    )
    .join('\n\n');

  vscode.window.showInformationMessage(
    `Parameter Mappings (${sourceFormat} -> ${targetFormat.label}):\n\n${mappingsText}`,
    'Close'
  );
}

/**
 * Migrate MD parameters
 */
async function migrateMDParameters(context: vscode.ExtensionContext): Promise<void> {
  vscode.window.showInformationMessage(
    'MD/Optimization Parameter Migration\n\n' +
    'This tool converts MD and optimization parameters:\n' +
    '- Time step (POTIM/dt/TIMESTEP)\n' +
    '- Temperature (TEBEG/TEMP/T)\n' +
    '- Pressure (PSTRESS/PRESS)\n' +
    '- Convergence criteria (EDIFF/conv_thr)\n\n' +
    'Note: MD parameter migration will be implemented in Phase 3.4',
    'OK'
  );
}

/**
 * Show migration validation report
 */
async function showMigrationValidation(
  context: vscode.ExtensionContext
): Promise<void> {
  vscode.window.showInformationMessage(
    'Migration Validation\n\n' +
    'Validation checks include:\n' +
    '- Atom count consistency\n' +
    '- Element composition\n' +
    '- Cell vector preservation\n' +
    '- Position accuracy\n\n' +
    'To validate a migration, run the migration with validation enabled.',
    'OK'
  );
}
