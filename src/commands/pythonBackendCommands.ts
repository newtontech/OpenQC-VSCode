/**
 * VS Code commands for the Scientific Python backend.
 *
 * Registers:
 * - openqc.checkPythonBackend
 * - openqc.configurePythonPath
 * - openqc.showPythonDiagnostics
 *
 * @module commands/pythonBackendCommands
 */

import * as vscode from 'vscode';
import { checkBackend, type BridgeResponse } from '../python/PythonBridge';
import type { BackendCheckResult } from '../python/PythonBackendStatus';
import { Logger } from '../utils/Logger';

const logger = Logger.getInstance();

// ---------------------------------------------------------------------------
// Output channel
// ---------------------------------------------------------------------------

let _channel: vscode.OutputChannel | undefined;

function getChannel(): vscode.OutputChannel {
  if (!_channel) {
    _channel = vscode.window.createOutputChannel('OpenQC Python Backend');
  }
  return _channel;
}

// ---------------------------------------------------------------------------
// Command registration
// ---------------------------------------------------------------------------

export function registerPythonBackendCommands(context: vscode.ExtensionContext): void {
  context.subscriptions.push(
    vscode.commands.registerCommand('openqc.checkPythonBackend', () => runBackendCheck()),
    vscode.commands.registerCommand('openqc.configurePythonPath', () => configurePythonPath()),
    vscode.commands.registerCommand('openqc.showPythonDiagnostics', () => showDiagnostics())
  );
}

// ---------------------------------------------------------------------------
// Implementation
// ---------------------------------------------------------------------------

async function runBackendCheck(): Promise<void> {
  await vscode.window.withProgress(
    {
      location: vscode.ProgressLocation.Notification,
      title: 'Checking Scientific Python Backend...',
      cancellable: true,
    },
    async (_progress, token) => {
      const result = await checkBackend({ cancelToken: token });
      displayBackendResult(result);
    }
  );
}

function displayBackendResult(result: BridgeResponse<BackendCheckResult>): void {
  const channel = getChannel();
  channel.clear();

  if (!result.success || !result.data) {
    const msg = result.error?.message ?? 'Unknown error';
    channel.appendLine(`❌ Backend check failed: ${msg}`);

    if (result.error?.hint) {
      channel.appendLine(`💡 Hint: ${result.error.hint}`);
    }

    channel.show(true);
    vscode.window.showErrorMessage(`Python backend check failed: ${msg}`);
    return;
  }

  const data = result.data;
  channel.appendLine('✅ Scientific Python Backend Check');
  channel.appendLine('═'.repeat(40));
  channel.appendLine('');

  // Python
  channel.appendLine(`🐍 Python: ${data.python.executable}`);
  channel.appendLine(`   Version: ${data.python.version}`);
  if (data.python.platform) {
    channel.appendLine(`   Platform: ${data.python.platform}`);
  }
  channel.appendLine('');

  // Packages
  channel.appendLine('📦 Scientific Packages:');
  const packages = data.packages ?? {};
  const pkgEntries = Object.entries(packages);
  const availableCount = pkgEntries.filter(([, s]) => s.available).length;

  for (const [name, status] of pkgEntries) {
    if (status.available) {
      channel.appendLine(`   ✅ ${name}: ${status.version ?? 'installed'}`);
    } else {
      channel.appendLine(`   ❌ ${name}: not installed`);
      if (status.installHint) {
        channel.appendLine(`      Install: ${status.installHint}`);
      }
    }
  }
  channel.appendLine(`\n   ${availableCount}/${pkgEntries.length} packages available`);
  channel.appendLine('');

  // External tools
  channel.appendLine('🔧 External Tools:');
  const tools = data.externalTools ?? {};
  for (const [name, status] of Object.entries(tools)) {
    if (status.available) {
      channel.appendLine(`   ✅ ${name}: ${status.path}`);
    } else {
      channel.appendLine(`   ❌ ${name}: not found`);
    }
  }
  channel.appendLine('');

  logger.info(`Backend check: ${availableCount}/${pkgEntries.length} packages available`);
  channel.show(true);

  vscode.window.showInformationMessage(
    `Python backend: ${data.python.version} — ${availableCount}/${pkgEntries.length} scientific packages available`
  );
}

async function configurePythonPath(): Promise<void> {
  const current = vscode.workspace.getConfiguration('openqc').get<string>('pythonPath', 'python3');

  const result = await vscode.window.showInputBox({
    prompt: 'Path to Python executable for OpenQC scientific backend',
    value: current,
    placeHolder: 'python3',
    validateInput: (value: string) => {
      if (!value.trim()) {
        return 'Path cannot be empty';
      }
      return undefined;
    },
  });

  if (result !== undefined) {
    await vscode.workspace
      .getConfiguration('openqc')
      .update('pythonPath', result, vscode.ConfigurationTarget.Global);
    vscode.window.showInformationMessage(`Python path set to: ${result}`);
  }
}

async function showDiagnostics(): Promise<void> {
  await runBackendCheck();
}
