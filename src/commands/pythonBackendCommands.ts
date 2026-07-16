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
  const backendStatus = data.status ?? 'degraded';
  const capabilitySummary = summarizeCapabilities(data.capabilities);
  channel.appendLine('✅ Scientific Python Backend Check');
  channel.appendLine('═'.repeat(40));
  channel.appendLine(`Status: ${backendStatus.toUpperCase()}`);
  if (data.statusDetail) {
    channel.appendLine(`Detail: ${data.statusDetail}`);
  }
  if (capabilitySummary) {
    channel.appendLine(`Capabilities: ${capabilitySummary}`);
  }
  if (data.missingPackages && data.missingPackages.length > 0) {
    channel.appendLine(`Missing core packages: ${data.missingPackages.join(', ')}`);
  }
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

  // Capabilities
  if (data.capabilities && Object.keys(data.capabilities).length > 0) {
    channel.appendLine('🧭 Feature Readiness:');
    for (const [name, capability] of Object.entries(data.capabilities)) {
      const icon = getCapabilityIcon(capability.status);
      const label = capability.label ?? humanizeCapabilityName(name);
      channel.appendLine(`   ${icon} ${label}: ${formatCapabilityStatus(capability.status)}`);
      channel.appendLine(`      ${capability.detail}`);
      if (capability.requires && capability.requires.length > 0) {
        channel.appendLine(`      Install for full support: ${capability.requires.join(', ')}`);
      }
    }
    channel.appendLine('');
  }

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

  const message = formatBackendNotification(data, availableCount, pkgEntries.length);
  if (backendStatus === 'installed') {
    vscode.window.showInformationMessage(message);
  } else if (backendStatus === 'degraded') {
    vscode.window.showWarningMessage(message);
  } else {
    vscode.window.showErrorMessage(message);
  }
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

function summarizeCapabilities(
  capabilities: BackendCheckResult['capabilities']
): string | undefined {
  if (!capabilities || Object.keys(capabilities).length === 0) {
    return undefined;
  }

  const counts = { available: 0, degraded: 0, missing: 0 };
  for (const capability of Object.values(capabilities)) {
    counts[capability.status] += 1;
  }

  return `${counts.available} available, ${counts.degraded} degraded, ${counts.missing} missing`;
}

function getCapabilityIcon(status: string): string {
  if (status === 'available') {
    return '✅';
  }
  if (status === 'degraded') {
    return '⚠️';
  }
  return '❌';
}

function formatCapabilityStatus(status: string): string {
  return status.toUpperCase();
}

function humanizeCapabilityName(name: string): string {
  return name.replace(/([a-z0-9])([A-Z])/g, '$1 $2').replace(/^./, first => first.toUpperCase());
}

function formatBackendNotification(
  data: BackendCheckResult,
  availableCount: number,
  totalCount: number
): string {
  const backendStatus = data.status ?? 'degraded';
  const packageSummary = `${availableCount}/${totalCount} scientific packages available`;

  if (backendStatus === 'installed') {
    return `Python backend installed: ${data.python.version} — ${packageSummary}`;
  }

  if (data.missingPackages && data.missingPackages.length > 0) {
    return `Python backend ${backendStatus}: missing ${data.missingPackages.join(
      ', '
    )}. ${packageSummary}`;
  }

  if (data.statusDetail) {
    return `Python backend ${backendStatus}: ${data.statusDetail}`;
  }

  return `Python backend ${backendStatus}: ${packageSummary}`;
}
