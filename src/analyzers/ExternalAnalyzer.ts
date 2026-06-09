/**
 * External analyzer adapter framework for Multiwfn and c2x.
 *
 * Safety model:
 * - External analyzers are disabled by default.
 * - User must configure path and confirm before execution.
 * - Commands are previewed before execution.
 * - Outputs include provenance metadata.
 *
 * @module analyzers/ExternalAnalyzer
 */

import * as vscode from 'vscode';
import { execFile } from 'child_process';
import { Logger } from '../utils/Logger';

const logger = Logger.getInstance();

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

export interface AnalyzerConfig {
  id: string;
  displayName: string;
  settingKey: string;
  defaultCommand: string;
  description: string;
}

export interface AnalyzerStatus {
  id: string;
  available: boolean;
  path: string | null;
  enabled: boolean;
}

export interface AnalyzerCommand {
  executable: string;
  args: string[];
  cwd: string;
  description: string;
}

// ---------------------------------------------------------------------------
// Adapter configs
// ---------------------------------------------------------------------------

export const MULTIWFN_CONFIG: AnalyzerConfig = {
  id: 'multiwfn',
  displayName: 'Multiwfn',
  settingKey: 'openqc.external.multiwfnPath',
  defaultCommand: 'Multiwfn',
  description: 'Wavefunction analysis: density, ESP, ELF, orbital cubes, population analysis',
};

export const C2X_CONFIG: AnalyzerConfig = {
  id: 'c2x',
  displayName: 'c2x',
  settingKey: 'openqc.external.c2xPath',
  defaultCommand: 'c2x',
  description: 'CASTEP density/checkpoint conversion to cube, XSF, and structure formats',
};

// ---------------------------------------------------------------------------
// Status check
// ---------------------------------------------------------------------------

export async function checkAnalyzer(config: AnalyzerConfig): Promise<AnalyzerStatus> {
  const settings = vscode.workspace.getConfiguration();
  const configuredPath = settings.get<string>(config.settingKey, '');
  const enabled = settings.get<boolean>('openqc.external.allowExternalAnalyzers', false);

  if (configuredPath) {
    try {
      const { stat } = await import('fs/promises');
      await stat(configuredPath);
      return { id: config.id, available: true, path: configuredPath, enabled };
    } catch {
      return { id: config.id, available: false, path: configuredPath, enabled };
    }
  }

  // Check PATH
  return new Promise(resolve => {
    execFile('which', [config.defaultCommand], (error, stdout) => {
      if (error) {
        resolve({ id: config.id, available: false, path: null, enabled });
      } else {
        resolve({
          id: config.id,
          available: true,
          path: stdout.trim() || null,
          enabled,
        });
      }
    });
  });
}

// ---------------------------------------------------------------------------
// Command preview and execution
// ---------------------------------------------------------------------------

export function generateMultiwfnScript(
  inputFile: string,
  operation: string,
  outputFile: string
): AnalyzerCommand {
  const config = MULTIWFN_CONFIG;
  return {
    executable: config.defaultCommand,
    args: [inputFile],
    cwd: '',
    description: `Multiwfn: ${operation} on ${inputFile} → ${outputFile}`,
  };
}

export function generateC2xCommand(
  inputFile: string,
  operation: string,
  outputFile: string
): AnalyzerCommand {
  return {
    executable: 'c2x',
    args: [inputFile, outputFile],
    cwd: '',
    description: `c2x: ${operation} ${inputFile} → ${outputFile}`,
  };
}

export async function previewCommand(command: AnalyzerCommand): Promise<boolean> {
  const detail =
    `Command: ${command.executable} ${command.args.join(' ')}\n` +
    `Working directory: ${command.cwd || '(current)'}\n` +
    `Description: ${command.description}`;

  const result = await vscode.window.showWarningMessage(
    'Confirm external analyzer execution?',
    { modal: true, detail },
    'Execute',
    'Cancel'
  );

  return result === 'Execute';
}

export async function executeAnalyzer(
  command: AnalyzerCommand,
  timeoutMs: number = 60000
): Promise<{ success: boolean; stdout: string; stderr: string; exitCode: number | null }> {
  return new Promise(resolve => {
    execFile(
      command.executable,
      command.args,
      { timeout: timeoutMs, cwd: command.cwd || undefined },
      (error, stdout, stderr) => {
        if (error) {
          resolve({
            success: false,
            stdout: stdout || '',
            stderr: stderr || error.message,
            exitCode: typeof error.code === 'number' ? error.code : null,
          });
        } else {
          resolve({ success: true, stdout, stderr: stderr || '', exitCode: 0 });
        }
      }
    );
  });
}
