/**
 * Migration Commands - VSCode Extension Commands
 *
 * Provides VSCode commands for structure migration, k-point migration,
 * parameter migration, and MD workflow migration.
 */

import * as path from 'path';
import * as vscode from 'vscode';
import { StructureMigration, quickMigrateStructure } from '../utils/migration/structure';
import { KpointMigration, quickMigrateKpoints } from '../utils/migration/kpoints';
import { ASEFormat } from '../ase/ASEConverter';
import { ParameterConverter } from '../utils/migration/parameterConverter';
import { MDWorkflowConverter, quickMigrateMDWorkflow } from '../utils/migration/mdWorkflow';

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

  // MD/Optimization workflow migration (Phase 3.4)
  const migrateMDWorkflowCommand = vscode.commands.registerCommand(
    'OpenQC.migrateMDWorkflow',
    async () => {
      await quickMigrateMDWorkflow(context);
    }
  );
  context.subscriptions.push(migrateMDWorkflowCommand);

  // Legacy MD parameters migration (for backward compatibility)
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
  const converter = new ParameterConverter();

  // Show progress
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Analyzing input file...',
      cancellable: false,
    },
    async () => {
      // Extract parameters
      const extraction = converter.extractParameters(sourcePath);

      if (!extraction.success || !extraction.data) {
        vscode.window.showErrorMessage(
          `Failed to extract parameters: ${extraction.error}`
        );
        return;
      }

      // Show extracted parameters
      const paramsList = Object.entries(extraction.data.parameters)
        .map(([key, value]) => `${key} = ${value}`)
        .join('\n');

      const action = await vscode.window.showInformationMessage(
        `Extracted ${Object.keys(extraction.data.parameters).length} parameters from ${extraction.data.format.toUpperCase()} file`,
        'Convert to Target Format',
        'View Parameters',
        'Copy to Clipboard'
      );

      if (action === 'Convert to Target Format') {
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

        // Convert parameters
        const conversion = await converter.convertFile(
          sourcePath,
          targetFormat.value
        );

        if (!conversion.success || !conversion.target) {
          vscode.window.showErrorMessage(
            `Conversion failed: ${conversion.error}`
          );
          return;
        }

        // Show converted parameters
        const convertedList = Object.entries(conversion.target.parameters)
          .map(([key, value]) => `${key} = ${value}`)
          .join('\n');

        const resultAction = await vscode.window.showInformationMessage(
          `Converted to ${targetFormat.label}:\n${convertedList}${conversion.target.unmapped.length > 0 ? `\n\nUnmapped: ${conversion.target.unmapped.join(', ')}` : ''}`,
          'Copy Converted Parameters',
          'View Warnings'
        );

        if (resultAction === 'Copy Converted Parameters') {
          vscode.env.clipboard.writeText(convertedList);
          vscode.window.showInformationMessage('Parameters copied to clipboard');
        } else if (resultAction === 'View Warnings' && conversion.warnings.length > 0) {
          vscode.window.showWarningMessage(conversion.warnings.join('\n'));
        }
      } else if (action === 'View Parameters') {
        // Create output channel to show parameters
        const channel = vscode.window.createOutputChannel('OpenQC Parameters');
        channel.clear();
        channel.appendLine(`Extracted Parameters (${extraction.data.format.toUpperCase()})`);
        channel.appendLine('='.repeat(50));
        channel.appendLine(paramsList);
        channel.appendLine('');
        channel.appendLine(`Total: ${Object.keys(extraction.data.parameters).length} parameters`);
        channel.show();
      } else if (action === 'Copy to Clipboard') {
        vscode.env.clipboard.writeText(paramsList);
        vscode.window.showInformationMessage('Parameters copied to clipboard');
      }
    }
  );
}

/**
 * Legacy MD parameters migration (for backward compatibility)
 * Use OpenQC.migrateMDWorkflow for comprehensive MD/Opt workflow migration
 */
async function migrateMDParameters(context: vscode.ExtensionContext): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const sourcePath = editor.document.uri.fsPath;
  const converter = new ParameterConverter();

  // Extract parameters
  const extraction = converter.extractParameters(sourcePath);

  if (!extraction.success || !extraction.data) {
    vscode.window.showErrorMessage(
      `Failed to extract parameters: ${extraction.error}`
    );
    return;
  }

  // Filter MD-related parameters
  const mdParams: Record<string, any> = {};
  const mdKeywords = ['POTIM', 'TIMESTEP', 'DT', 'TEBEG', 'TEMP', 'T', 'PSTRESS', 'PRESS', 'NSW', 'MD'];

  Object.entries(extraction.data.parameters).forEach(([key, value]) => {
    if (mdKeywords.some(kw => key.toUpperCase().includes(kw))) {
      mdParams[key] = value;
    }
  });

  if (Object.keys(mdParams).length === 0) {
    vscode.window.showInformationMessage(
      'No MD parameters found in this file'
    );
    return;
  }

  const mdList = Object.entries(mdParams)
    .map(([key, value]) => `${key} = ${value}`)
    .join('\n');

  const action = await vscode.window.showInformationMessage(
    `Found ${Object.keys(mdParams).length} MD parameters:\n${mdList}`,
    'Convert to Target Format',
    'Copy to Clipboard',
    'Use MD Workflow Migration'
  );

  if (action === 'Use MD Workflow Migration') {
    // Suggest using the comprehensive MD workflow migration
    vscode.window.showInformationMessage(
      'MD Workflow Migration provides:\n' +
      '- Ensemble type detection (NVE/NVT/NPT/NPH)\n' +
      '- Thermostat and barostat conversion\n' +
      '- Optimization algorithm mapping\n' +
      '- Comprehensive parameter generation\n\n' +
      'Use: OpenQC: Migrate MD Workflow command'
    );
  } else if (action === 'Convert to Target Format') {
    // Select target format
    const targetFormat = await vscode.window.showQuickPick(
      ['vasp', 'cp2k', 'qe'].map(code => ({
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

    // Convert parameters
    const conversion = await converter.convertFile(
      sourcePath,
      targetFormat.value
    );

    if (!conversion.success || !conversion.target) {
      vscode.window.showErrorMessage(
        `Conversion failed: ${conversion.error}`
      );
      return;
    }

    // Show converted MD parameters
    const convertedMD: Record<string, any> = {};
    Object.entries(conversion.target.parameters).forEach(([key, value]) => {
      if (mdKeywords.some(kw => key.toUpperCase().includes(kw))) {
        convertedMD[key] = value;
      }
    });

    const convertedList = Object.entries(convertedMD)
      .map(([key, value]) => `${key} = ${value}`)
      .join('\n');

    vscode.window.showInformationMessage(
      `Converted MD parameters to ${targetFormat.label}:\n${convertedList}`
    );
  } else if (action === 'Copy to Clipboard') {
    vscode.env.clipboard.writeText(mdList);
    vscode.window.showInformationMessage('MD parameters copied to clipboard');
  }
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
    'To validate a migration, run the migration with validation enabled.\n\n' +
    'For MD/Opt workflows:\n' +
    '- Ensemble type correctness\n' +
    '- Thermostat/barostat parameter mapping\n' +
    '- Convergence criteria units',
    'OK'
  );
}
