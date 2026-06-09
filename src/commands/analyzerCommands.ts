/**
 * VS Code commands for external analyzer adapters.
 *
 * @module commands/analyzerCommands
 */

import * as vscode from 'vscode';
import {
  checkAnalyzer,
  MULTIWFN_CONFIG,
  C2X_CONFIG,
  generateMultiwfnScript,
  generateC2xCommand,
  previewCommand,
  executeAnalyzer,
} from '../analyzers/ExternalAnalyzer';

export function registerAnalyzerCommands(context: vscode.ExtensionContext): void {
  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.checkExternalAnalyzers', () => checkAnalyzers()),
    vscode.commands.registerCommand('openqc.generateMultiwfnScript', () =>
      generateMultiwfnScriptCommand()
    ),
    vscode.commands.registerCommand('openqc.runMultiwfnAnalysis', () => runMultiwfnCommand()),
    vscode.commands.registerCommand('openqc.convertDensityC2x', () => convertWithC2xCommand())
  );
}

async function checkAnalyzers(): Promise<void> {
  const [multiwfn, c2x] = await Promise.all([
    checkAnalyzer(MULTIWFN_CONFIG),
    checkAnalyzer(C2X_CONFIG),
  ]);

  const channel = vscode.window.createOutputChannel('OpenQC External Analyzers');
  channel.clear();
  channel.appendLine('External Analyzer Status');
  channel.appendLine('═'.repeat(40));
  channel.appendLine('');

  for (const status of [multiwfn, c2x]) {
    const icon = status.available ? '✅' : '❌';
    const path = status.path ?? 'not found';
    const enabled = status.enabled ? 'enabled' : 'disabled';
    channel.appendLine(`${icon} ${status.id}: ${path} (${enabled})`);
  }

  const enabled = vscode.workspace
    .getConfiguration('openqc.external')
    .get<boolean>('allowExternalAnalyzers', false);
  channel.appendLine(
    `\nAnalyzers are ${enabled ? 'ENABLED' : 'DISABLED'} (openqc.external.allowExternalAnalyzers)`
  );
  channel.show(true);
}

async function generateMultiwfnScriptCommand(): Promise<void> {
  const inputFile = await vscode.window.showInputBox({
    prompt: 'Input file path (.fchk, .wfn, .molden, .cube)',
    placeHolder: '/path/to/file.fchk',
  });
  if (!inputFile) return;

  const operation = await vscode.window.showQuickPick(
    [
      'Electron density cube',
      'Electrostatic potential cube',
      'ELF cube',
      'Population analysis',
      'Orbital cube',
    ],
    { placeHolder: 'Select analysis type' }
  );
  if (!operation) return;

  const outputFile = `${inputFile}.${operation.toLowerCase().replace(/\s+/g, '_')}.out`;
  const command = generateMultiwfnScript(inputFile, operation, outputFile);

  const doc = await vscode.workspace.openTextDocument({
    content: `# Multiwfn Analysis Script\n# ${command.description}\n\n${command.executable} ${command.args.join(' ')}\n`,
    language: 'shellscript',
  });
  await vscode.window.showTextDocument(doc);

  vscode.window.showInformationMessage(
    `Script generated. Run manually after review. Analyzers are disabled by default for safety.`
  );
}

async function runMultiwfnCommand(): Promise<void> {
  vscode.window.showWarningMessage(
    'External analyzers must be enabled in settings before execution. Set openqc.external.allowExternalAnalyzers to true.'
  );
}

async function convertWithC2xCommand(): Promise<void> {
  vscode.window.showInformationMessage(
    'c2x conversion requires CASTEP checkpoint or density files. Configure path in openqc.external.c2xPath.'
  );
}
