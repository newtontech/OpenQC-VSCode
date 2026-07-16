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
  OPENBABEL_CONFIG,
  generateMultiwfnScript,
  generateC2xCommand,
  generateOpenBabelCommand,
  previewCommand,
  executeAnalyzer,
} from '../analyzers/ExternalAnalyzer';
import type { AnalyzerCommand } from '../analyzers/ExternalAnalyzer';

export function registerAnalyzerCommands(context: vscode.ExtensionContext): void {
  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.checkExternalAnalyzers', () => checkAnalyzers()),
    vscode.commands.registerCommand('openqc.generateMultiwfnScript', () =>
      generateMultiwfnScriptCommand()
    ),
    vscode.commands.registerCommand('openqc.runMultiwfnAnalysis', () => runMultiwfnCommand()),
    vscode.commands.registerCommand('openqc.convertDensityC2x', () => convertWithC2xCommand()),
    vscode.commands.registerCommand('openqc.convertStructureOpenBabel', () =>
      convertWithOpenBabelCommand()
    )
  );
}

async function checkAnalyzers(): Promise<void> {
  const [multiwfn, c2x, openbabel] = await Promise.all([
    checkAnalyzer(MULTIWFN_CONFIG),
    checkAnalyzer(C2X_CONFIG),
    checkAnalyzer(OPENBABEL_CONFIG),
  ]);

  const channel = vscode.window.createOutputChannel('OpenQC External Analyzers');
  channel.clear();
  channel.appendLine('External Analyzer Status');
  channel.appendLine('═'.repeat(40));
  channel.appendLine('');

  for (const status of [multiwfn, c2x, openbabel]) {
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
  if (!inputFile) {
    return;
  }

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
  if (!operation) {
    return;
  }

  const outputFile = defaultMultiwfnOutputPath(inputFile, operation);
  const command = generateMultiwfnScript(inputFile, operation, outputFile);

  const doc = await vscode.workspace.openTextDocument({
    content: renderShellCommand(command),
    language: 'shellscript',
  });
  await vscode.window.showTextDocument(doc);

  vscode.window.showInformationMessage(
    `Script generated. Run manually after review. Analyzers are disabled by default for safety.`
  );
}

async function runMultiwfnCommand(): Promise<void> {
  const status = await checkAnalyzer(MULTIWFN_CONFIG);
  if (!ensureAnalyzerReady(status, MULTIWFN_CONFIG.displayName)) {
    return;
  }

  const inputFile = await vscode.window.showInputBox({
    prompt: 'Input file path for Multiwfn (.fchk, .wfn, .molden, .cube)',
    placeHolder: '/path/to/file.fchk',
  });
  if (!inputFile) {
    return;
  }

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
  if (!operation) {
    return;
  }

  const outputFile = defaultMultiwfnOutputPath(inputFile, operation);
  const command = generateMultiwfnScript(inputFile, operation, outputFile);
  command.executable = status.path ?? command.executable;

  await runAnalyzerCommand(command);
}

function defaultMultiwfnOutputPath(inputFile: string, operation: string): string {
  const suffix = operation.toLowerCase().replace(/\s+/g, '_');
  const extension = operation.toLowerCase().includes('cube') ? 'cube' : 'out';
  return `${inputFile}.${suffix}.${extension}`;
}

async function convertWithC2xCommand(): Promise<void> {
  const status = await checkAnalyzer(C2X_CONFIG);
  if (!ensureAnalyzerReady(status, C2X_CONFIG.displayName)) {
    return;
  }

  const inputFile = await vscode.window.showInputBox({
    prompt: 'Input density/checkpoint file for c2x',
    placeHolder: '/path/to/file.check',
  });
  if (!inputFile) {
    return;
  }

  const outputFile = await vscode.window.showInputBox({
    prompt: 'Output file path for c2x conversion',
    value: `${inputFile}.cube`,
  });
  if (!outputFile) {
    return;
  }

  const command = generateC2xCommand(inputFile, 'convert', outputFile);
  command.executable = status.path ?? command.executable;

  await runAnalyzerCommand(command);
}

async function convertWithOpenBabelCommand(): Promise<void> {
  const status = await checkAnalyzer(OPENBABEL_CONFIG);
  if (!ensureAnalyzerReady(status, OPENBABEL_CONFIG.displayName)) {
    return;
  }

  const inputFile = await vscode.window.showInputBox({
    prompt: 'Input structure or molecule file for Open Babel',
    placeHolder: '/path/to/structure.xyz',
  });
  if (!inputFile) {
    return;
  }

  const outputFile = await vscode.window.showInputBox({
    prompt: 'Output structure or molecule file for Open Babel',
    value: `${inputFile}.pdb`,
  });
  if (!outputFile) {
    return;
  }

  const command = generateOpenBabelCommand(inputFile, outputFile);
  command.executable = status.path ?? command.executable;

  await runAnalyzerCommand(command);
}

function ensureAnalyzerReady(
  status: { available: boolean; enabled: boolean; path: string | null },
  displayName: string
): boolean {
  if (!status.enabled) {
    vscode.window.showWarningMessage(
      `External analyzers are disabled. Enable openqc.external.allowExternalAnalyzers before running ${displayName}.`
    );
    return false;
  }

  if (!status.available || !status.path) {
    vscode.window.showErrorMessage(
      `${displayName} executable was not found. Configure the corresponding openqc.external.*Path setting or add it to PATH.`
    );
    return false;
  }

  return true;
}

async function runAnalyzerCommand(command: AnalyzerCommand): Promise<void> {
  if (!(await previewCommand(command))) {
    return;
  }

  const timeoutMs = vscode.workspace
    .getConfiguration('openqc.external')
    .get<number>('timeoutMs', 60000);

  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: command.description,
      cancellable: false,
    },
    async () => {
      const result = await executeAnalyzer(command, timeoutMs);
      const channel = vscode.window.createOutputChannel('OpenQC External Analyzer Run');
      channel.clear();
      channel.appendLine(command.description);
      channel.appendLine(`Command: ${command.executable} ${command.args.join(' ')}`);
      if (command.expectedOutputPath) {
        channel.appendLine(`Expected output: ${command.expectedOutputPath}`);
      }
      channel.appendLine(`Exit code: ${result.exitCode}`);
      channel.appendLine('');
      channel.appendLine('--- stdout ---');
      channel.appendLine(result.stdout || '(empty)');
      channel.appendLine('--- stderr ---');
      channel.appendLine(result.stderr || '(empty)');
      channel.show(true);

      if (result.success) {
        vscode.window.showInformationMessage('External analyzer completed successfully');
      } else {
        vscode.window.showErrorMessage(
          'External analyzer failed; see OpenQC External Analyzer Run output'
        );
      }
    }
  );
}

function renderShellCommand(command: AnalyzerCommand): string {
  const commandLine = `${command.executable} ${command.args.join(' ')}`;
  if (!command.stdin) {
    return `# ${command.description}\n\n${commandLine}\n`;
  }

  return `# ${command.description}\n\n${commandLine} <<'EOF'\n${command.stdin}EOF\n`;
}
