/**
 * ASE Integration Commands
 *
 * VSCode commands for ASE format conversion and migration.
 */

import * as vscode from 'vscode';
import { ASEConverter, ASEFormat, ConversionResult } from './ASEConverter';

/**
 * Register ASE commands
 */
export function registerASECommands(context: vscode.ExtensionContext): void {
  const converter = new ASEConverter(context);

  // Command: Convert to ASE Atoms
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.convertToASE',
      async () => await convertToASE(converter)
    )
  );

  // Command: Convert from ASE Atoms
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.convertFromASE',
      async () => await convertFromASE(converter)
    )
  );

  // Command: Migrate to target format
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.migrateFormat',
      async () => await migrateFormat(converter)
    )
  );

  // Command: Quick convert (VASP -> CP2K, etc.)
  context.subscriptions.push(
    vscode.commands.registerCommand(
      'openqc.quickConvert',
      async () => await quickConvert(converter)
    )
  );
}

/**
 * Convert current file to ASE Atoms
 */
async function convertToASE(converter: ASEConverter): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const document = editor.document;
  const filepath = document.uri.fsPath;

  // Show progress
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Converting to ASE Atoms',
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Reading file...' });

      const result = await converter.readToAtoms(filepath);

      if (!result.success) {
        vscode.window.showErrorMessage(`Conversion failed: ${result.error}`);
        return;
      }

      progress.report({ message: 'Displaying results...' });

      // Display results in a new document
      const atomsInfo = formatAtomsInfo(result);
      const doc = await vscode.workspace.openTextDocument({
        content: atomsInfo,
        language: 'json',
      });
      await vscode.window.showTextDocument(doc);

      // Show success message
      const formula = result.atoms?.chemical_symbols
        ? getChemicalFormula(result.atoms.chemical_symbols)
        : 'Unknown';
      vscode.window.showInformationMessage(
        `Successfully converted to ASE Atoms: ${formula} (${result.metadata.natoms} atoms)`
      );
    }
  );
}

/**
 * Convert ASE Atoms to target format
 */
async function convertFromASE(converter: ASEConverter): Promise<void> {
  // This would typically be called after convertToASE
  // For now, show a quick pick of available formats
  const formats = converter.getSupportedFormats();
  const items = Object.entries(formats).map(([key, value]) => ({
    label: value.name,
    description: value.description,
    detail: `Extensions: ${value.extensions.join(', ')}`,
    format: key as ASEFormat,
  }));

  const selected = await vscode.window.showQuickPick(items, {
    placeHolder: 'Select target format',
  });

  if (!selected) {
    return;
  }

  // Ask for output file path
  const outputPath = await vscode.window.showInputBox({
    prompt: 'Enter output file path',
    value: `output${formats[selected.format].extensions[0]}`,
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a file path';
      }
      return null;
    },
  });

  if (!outputPath) {
    return;
  }

  vscode.window.showInformationMessage(
    `Converting to ${selected.format}... (Feature coming soon - need atoms data)`
  );
}

/**
 * Migrate from one format to another
 */
async function migrateFormat(converter: ASEConverter): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const inputPath = editor.document.uri.fsPath;

  // Select target format
  const formats = converter.getSupportedFormats();
  const items = Object.entries(formats).map(([key, value]) => ({
    label: value.name,
    description: value.description,
    detail: `Extensions: ${value.extensions.join(', ')}`,
    format: key as ASEFormat,
  }));

  const selected = await vscode.window.showQuickPick(items, {
    placeHolder: 'Select target format for migration',
  });

  if (!selected) {
    return;
  }

  // Suggest output filename
  const inputBasename = inputPath.replace(/\.[^/.]+$/, '');
  const suggestedOutput = `${inputBasename}_${selected.format}${
    formats[selected.format].extensions[0]
  }`;

  const outputPath = await vscode.window.showInputBox({
    prompt: 'Enter output file path',
    value: suggestedOutput,
    validateInput: value => {
      if (!value || value.trim().length === 0) {
        return 'Please enter a file path';
      }
      return null;
    },
  });

  if (!outputPath) {
    return;
  }

  // Perform conversion
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Migrating to ${selected.format}`,
      cancellable: false,
    },
    async progress => {
      progress.report({ message: 'Converting format...' });

      const result = await converter.convertFormat(
        inputPath,
        outputPath,
        undefined,
        selected.format
      );

      if (!result.success) {
        vscode.window.showErrorMessage(`Migration failed: ${result.error}`);
        return;
      }

      progress.report({ message: 'Opening result...' });

      // Open the converted file
      const doc = await vscode.workspace.openTextDocument(outputPath);
      await vscode.window.showTextDocument(doc);

      vscode.window.showInformationMessage(
        `Successfully migrated to ${selected.format}: ${outputPath}`
      );
    }
  );
}

/**
 * Quick convert with common presets
 */
async function quickConvert(converter: ASEConverter): Promise<void> {
  const presets = [
    { label: 'VASP → CP2K', input: ASEFormat.VASP, output: ASEFormat.CP2K },
    { label: 'VASP → Quantum ESPRESSO', input: ASEFormat.VASP, output: ASEFormat.QE },
    { label: 'CP2K → VASP', input: ASEFormat.CP2K, output: ASEFormat.VASP },
    { label: 'Gaussian → ORCA', input: ASEFormat.Gaussian, output: ASEFormat.ORCA },
    { label: 'ORCA → Gaussian', input: ASEFormat.ORCA, output: ASEFormat.Gaussian },
    { label: 'XYZ → VASP', input: ASEFormat.XYZ, output: ASEFormat.VASP },
    { label: 'CIF → VASP', input: ASEFormat.CIF, output: ASEFormat.VASP },
  ];

  const selected = await vscode.window.showQuickPick(
    presets.map(p => p),
    { placeHolder: 'Select conversion preset' }
  );

  if (!selected) {
    return;
  }

  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active editor');
    return;
  }

  const inputPath = editor.document.uri.fsPath;
  const inputBasename = inputPath.replace(/\.[^/.]+$/, '');
  const outputPath = `${inputBasename}_${selected.output}${getExtension(selected.output)}`;

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Converting: ${selected.label}`,
      cancellable: false,
    },
    async progress => {
      const result = await converter.convertFormat(
        inputPath,
        outputPath,
        selected.input,
        selected.output
      );

      if (!result.success) {
        vscode.window.showErrorMessage(`Conversion failed: ${result.error}`);
        return;
      }

      const doc = await vscode.workspace.openTextDocument(outputPath);
      await vscode.window.showTextDocument(doc);

      vscode.window.showInformationMessage(`Successfully converted: ${outputPath}`);
    }
  );
}

/**
 * Format ASE Atoms information for display
 */
function formatAtomsInfo(result: ConversionResult): string {
  if (!result.atoms) {
    return JSON.stringify(result, null, 2);
  }

  const atoms = result.atoms;
  const formula = getChemicalFormula(atoms.chemical_symbols);

  const info = {
    summary: {
      formula: formula,
      natoms: atoms.chemical_symbols.length,
      periodic: atoms.pbc.some(p => p),
      elements: [...new Set(atoms.chemical_symbols)],
    },
    structure: {
      positions: atoms.positions,
      cell: atoms.cell,
      pbc: atoms.pbc,
    },
    metadata: result.metadata,
    warnings: result.warnings,
  };

  return JSON.stringify(info, null, 2);
}

/**
 * Get chemical formula from element symbols
 */
function getChemicalFormula(symbols: string[]): string {
  const counts: Record<string, number> = {};
  for (const symbol of symbols) {
    counts[symbol] = (counts[symbol] || 0) + 1;
  }

  return Object.entries(counts)
    .map(([element, count]) => `${element}${count > 1 ? count : ''}`)
    .join('');
}

/**
 * Get default extension for format
 */
function getExtension(format: ASEFormat): string {
  const extensions: Record<ASEFormat, string> = {
    [ASEFormat.VASP]: '.POSCAR',
    [ASEFormat.CP2K]: '.inp',
    [ASEFormat.QE]: '.in',
    [ASEFormat.Gaussian]: '.com',
    [ASEFormat.ORCA]: '.inp',
    [ASEFormat.NWChem]: '.nw',
    [ASEFormat.GAMESS]: '.inp',
    [ASEFormat.LAMMPS]: '.lmp',
    [ASEFormat.XYZ]: '.xyz',
    [ASEFormat.PDB]: '.pdb',
    [ASEFormat.CIF]: '.cif',
  };
  return extensions[format] || '.txt';
}
