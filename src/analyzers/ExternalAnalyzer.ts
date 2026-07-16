/**
 * External analyzer adapter framework for Multiwfn, c2x, and Open Babel.
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
import { constants } from 'fs';
import { access } from 'fs/promises';
import * as path from 'path';
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
  stdin?: string;
  expectedOutputPath?: string;
}

export interface PathLookupCommand {
  executable: string;
  args: string[];
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

export const OPENBABEL_CONFIG: AnalyzerConfig = {
  id: 'openbabel',
  displayName: 'Open Babel',
  settingKey: 'openqc.external.openBabelPath',
  defaultCommand: 'obabel',
  description: 'Structure and molecule format conversion via the Open Babel CLI',
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
      await access(configuredPath, constants.X_OK);
      return { id: config.id, available: true, path: configuredPath, enabled };
    } catch {
      return { id: config.id, available: false, path: configuredPath, enabled };
    }
  }

  // Check PATH
  const lookup = getPathLookupCommand(config.defaultCommand);
  return new Promise(resolve => {
    execFile(lookup.executable, lookup.args, (error, stdout) => {
      if (error) {
        resolve({ id: config.id, available: false, path: null, enabled });
      } else {
        const resolvedPath = firstResolvedPath(stdout);
        resolve({
          id: config.id,
          available: Boolean(resolvedPath),
          path: resolvedPath,
          enabled,
        });
      }
    });
  });
}

export function getPathLookupCommand(
  defaultCommand: string,
  platform: NodeJS.Platform = process.platform
): PathLookupCommand {
  if (platform === 'win32') {
    return { executable: 'where', args: [defaultCommand] };
  }
  return { executable: 'which', args: [defaultCommand] };
}

function firstResolvedPath(stdout: string): string | null {
  return (
    stdout
      .split(/\r?\n/)
      .map(line => line.trim())
      .filter(Boolean)[0] ?? null
  );
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
    stdin: generateMultiwfnStdin(operation, outputFile),
    expectedOutputPath: outputFile,
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
    expectedOutputPath: outputFile,
    description: `c2x: ${operation} ${inputFile} → ${outputFile}`,
  };
}

export function generateOpenBabelCommand(inputFile: string, outputFile: string): AnalyzerCommand {
  return {
    executable: 'obabel',
    args: [inputFile, '-O', outputFile],
    cwd: '',
    expectedOutputPath: outputFile,
    description: `Open Babel: convert ${inputFile} → ${outputFile}`,
  };
}

export async function previewCommand(command: AnalyzerCommand): Promise<boolean> {
  const detail =
    `Command: ${command.executable} ${command.args.join(' ')}\n` +
    `Working directory: ${command.cwd || '(current)'}\n` +
    (command.stdin ? `Standard input:\n${command.stdin}\n` : '') +
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
    const child = execFile(
      command.executable,
      command.args,
      { timeout: timeoutMs, cwd: command.cwd || undefined },
      async (error, stdout, stderr) => {
        if (error) {
          resolve({
            success: false,
            stdout: stdout || '',
            stderr: stderr || error.message,
            exitCode: typeof error.code === 'number' ? error.code : null,
          });
        } else {
          if (command.expectedOutputPath) {
            const expectedOutputPath = path.isAbsolute(command.expectedOutputPath)
              ? command.expectedOutputPath
              : path.join(command.cwd || process.cwd(), command.expectedOutputPath);
            try {
              await access(expectedOutputPath);
            } catch {
              resolve({
                success: false,
                stdout,
                stderr: stderr || `Expected analyzer output was not created: ${expectedOutputPath}`,
                exitCode: 0,
              });
              return;
            }
          }
          resolve({ success: true, stdout, stderr: stderr || '', exitCode: 0 });
        }
      }
    );
    if (command.stdin && child.stdin) {
      child.stdin.write(command.stdin);
      child.stdin.end();
    }
  });
}

function generateMultiwfnStdin(operation: string, outputFile: string): string {
  const normalized = operation.toLowerCase();
  if (normalized.includes('population')) {
    return ['7', '1', '0', 'q'].join('\n') + '\n';
  }

  if (normalized.includes('orbital')) {
    return ['200', '3', outputFile, '0', 'q'].join('\n') + '\n';
  }

  if (normalized.includes('electrostatic')) {
    return ['5', '12', '2', outputFile, '0', 'q'].join('\n') + '\n';
  }

  if (normalized.includes('elf')) {
    return ['5', '9', '2', outputFile, '0', 'q'].join('\n') + '\n';
  }

  return ['5', '1', '2', outputFile, '0', 'q'].join('\n') + '\n';
}
