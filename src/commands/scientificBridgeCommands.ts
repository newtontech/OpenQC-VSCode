/**
 * VS Code commands for scientific Python bridges.
 *
 * Registers commands for structure parsing, output analysis,
 * and software coverage expansion.
 *
 * @module commands/scientificBridgeCommands
 */

import * as vscode from 'vscode';
import { parseStructure, generateSupercell } from '../python/StructureBridge';
import { parseOutput } from '../results/OutputParserBridge';
import { Logger } from '../utils/Logger';

const logger = Logger.getInstance();

export function registerScientificBridgeCommands(context: vscode.ExtensionContext): void {
  context.subscriptions.push(
    // Structure bridge commands (#75)
    vscode.commands.registerCommand('openqc.parseStructurePython', () => parseStructureCommand()),
    vscode.commands.registerCommand('openqc.generateSupercell', () => generateSupercellCommand()),

    // Output parser commands (#76)
    vscode.commands.registerCommand('openqc.parseCalculationOutput', () => parseOutputCommand()),
    vscode.commands.registerCommand('openqc.showSCFConvergence', () => showSCFConvergenceCommand())
  );
}

// ---------------------------------------------------------------------------
// Structure commands
// ---------------------------------------------------------------------------

async function parseStructureCommand(): Promise<void> {
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active file');
    return;
  }

  const filePath = editor.document.uri.fsPath;
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Parsing structure with Python backend...',
      cancellable: true,
    },
    async (_progress, token) => {
      const result = await parseStructure(filePath, undefined, undefined);

      if (!result.success || !result.data) {
        const msg = result.error?.message ?? 'Unknown error';
        vscode.window.showErrorMessage(`Structure parse failed: ${msg}`);
        return;
      }

      const structure = result.data;
      const atomCount = structure.atoms?.length ?? 0;
      const kind = structure.kind ?? 'unknown';

      vscode.window.showInformationMessage(`Parsed ${atomCount} atoms (${kind}) from ${filePath}`);

      // Open structure JSON in a new editor
      const doc = await vscode.workspace.openTextDocument({
        content: JSON.stringify(structure, null, 2),
        language: 'json',
      });
      await vscode.window.showTextDocument(doc, vscode.ViewColumn.Beside);
    }
  );
}

async function generateSupercellCommand(): Promise<void> {
  // Get active file or ask user
  const editor = vscode.window.activeTextEditor;
  if (!editor) {
    vscode.window.showErrorMessage('No active file. Open a periodic structure file first.');
    return;
  }

  const result = await vscode.window.showInputBox({
    prompt: 'Supercell dimensions (e.g., 2 2 2)',
    placeHolder: '2 2 2',
    value: '2 2 2',
  });

  if (!result) return;

  const parts = result.trim().split(/\s+/).map(Number);
  if (parts.length !== 3 || parts.some(isNaN)) {
    vscode.window.showErrorMessage('Invalid dimensions. Use format: Na Nb Nc (e.g., 2 2 2)');
    return;
  }

  const [na, nb, nc] = parts;
  const maxAtoms = 10000;

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: `Generating ${na}x${nb}x${nc} supercell...`,
    },
    async () => {
      // First parse the structure
      const filePath = editor.document.uri.fsPath;
      const parseResult = await parseStructure(filePath, undefined, editor.document.getText());

      if (!parseResult.success || !parseResult.data) {
        vscode.window.showErrorMessage(`Parse failed: ${parseResult.error?.message}`);
        return;
      }

      const structure = parseResult.data;
      if (!structure.cell) {
        vscode.window.showErrorMessage('Structure is not periodic — cannot generate supercell');
        return;
      }

      const atomCount = structure.atoms.length * na * nb * nc;
      if (atomCount > maxAtoms) {
        vscode.window.showWarningMessage(
          `Supercell would have ${atomCount} atoms (max ${maxAtoms}). Reduce dimensions.`
        );
        return;
      }

      const scResult = await generateSupercell(structure, na, nb, nc);

      if (!scResult.success || !scResult.data) {
        vscode.window.showErrorMessage(`Supercell failed: ${scResult.error?.message}`);
        return;
      }

      vscode.window.showInformationMessage(
        `Supercell ${na}x${nb}x${nc}: ${scResult.data.atoms?.length ?? 0} atoms`
      );

      // Show in viewer
      const { OpenQCViewerPanel } = await import('../visualizers/OpenQCViewerPanel');
      OpenQCViewerPanel.createOrShow(
        vscode.Uri.parse(''),
        scResult.data,
        `Supercell ${na}x${nb}x${nc}`
      );
    }
  );
}

// ---------------------------------------------------------------------------
// Output parser commands
// ---------------------------------------------------------------------------

async function parseOutputCommand(): Promise<void> {
  const uris = await vscode.window.showOpenDialog({
    canSelectMany: false,
    title: 'Select calculation output file',
    filters: {
      'Output files': ['log', 'out', 'dat'],
      'All files': ['*'],
    },
  });

  if (!uris || uris.length === 0) return;

  const filePath = uris[0].fsPath;

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Parsing calculation output...',
    },
    async () => {
      const result = await parseOutput(filePath);

      if (!result.success || !result.data) {
        vscode.window.showErrorMessage(`Output parse failed: ${result.error?.message}`);
        return;
      }

      const results = result.data;
      const energy = results.finalEnergy
        ? `${results.finalEnergy.value.toFixed(4)} ${results.finalEnergy.unit}`
        : 'N/A';

      vscode.window.showInformationMessage(
        `Parsed output: ${results.software ?? 'unknown'}, energy=${energy}, SCF steps=${results.scfEnergies?.length ?? 0}`
      );

      // Show results JSON
      const doc = await vscode.workspace.openTextDocument({
        content: JSON.stringify(results, null, 2),
        language: 'json',
      });
      await vscode.window.showTextDocument(doc, vscode.ViewColumn.Beside);
    }
  );
}

async function showSCFConvergenceCommand(): Promise<void> {
  vscode.window.showInformationMessage(
    'Open a calculation output file and run "OpenQC: Parse Calculation Output" to view SCF convergence.'
  );
}
