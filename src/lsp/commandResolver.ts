/**
 * Cross-platform LSP Command Resolver
 *
 * Resolves language server startup commands without shell interpolation.
 * Supports PATH lookups, absolute paths, paths with spaces, and
 * Python-module-based servers.
 *
 * @module lsp/commandResolver
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/97
 */

import { execFile } from 'child_process';
import * as path from 'path';
import * as vscode from 'vscode';
import { LSPServerRegistryEntry } from './types';

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

/**
 * A fully resolved command that can be passed directly to `LanguageClient`.
 */
export type ResolvedLspCommand =
  | { kind: 'pathOrCommand'; command: string; args: string[]; env?: Record<string, string> }
  | {
      kind: 'pythonModule';
      python: string;
      module: string;
      args: string[];
      env?: Record<string, string>;
    };

/**
 * User-configurable overrides read from VS Code settings.
 */
export interface LspCommandOverrides {
  /** Deprecated alias: absolute path or executable name. */
  path?: string;
  /** Explicit command (executable name or absolute path). */
  command?: string;
  /** Arguments forwarded to the server process. */
  args?: string[];
  /** Extra environment variables for the server process. */
  env?: Record<string, string>;
}

// ---------------------------------------------------------------------------
// Resolver
// ---------------------------------------------------------------------------

/**
 * Resolve an LSP startup command by merging registry defaults with
 * user-provided config overrides.
 *
 * Resolution order (last wins):
 *   1. Registry `executable` (default command name)
 *   2. Deprecated `openqc.lsp.<id>.path` setting
 *   3. New `openqc.lsp.<id>.command` setting
 *
 * Args come from:
 *   1. Registry default `['--stdio']` if the registry entry provides no args
 *   2. User `openqc.lsp.<id>.args` override
 *
 * @param entry     - Bundled registry entry (source of defaults).
 * @param overrides - User settings read from VS Code configuration.
 */
export function resolveLspCommand(
  entry: LSPServerRegistryEntry,
  overrides: LspCommandOverrides
): ResolvedLspCommand {
  const command = overrides.command || overrides.path || entry.executable;
  const args = overrides.args ?? [...(entry.args ?? ['--stdio'])];
  const env = overrides.env;

  return {
    kind: 'pathOrCommand',
    command,
    args,
    env: Object.keys(env || {}).length > 0 ? env : undefined,
  };
}

/**
 * Read user-facing LSP command overrides for the given language ID from
 * VS Code configuration.
 *
 * Looks up:
 *   - `openqc.lsp.<languageId>.path`   (deprecated)
 *   - `openqc.lsp.<languageId>.command`
 *   - `openqc.lsp.<languageId>.args`
 *   - `openqc.lsp.<languageId>.env`
 */
export function readCommandOverrides(
  config: vscode.WorkspaceConfiguration,
  languageId: string
): LspCommandOverrides {
  return {
    path: config.get<string | undefined>(`${languageId}.path`, undefined),
    command: config.get<string | undefined>(`${languageId}.command`, undefined),
    args: config.get<string[] | undefined>(`${languageId}.args`, undefined),
    env: config.get<Record<string, string> | undefined>(`${languageId}.env`, undefined),
  };
}

// ---------------------------------------------------------------------------
// Executable existence check (cross-platform)
// ---------------------------------------------------------------------------

/**
 * Check whether an executable is resolvable on the current system.
 *
 * - On Windows, uses `where <command>`.
 * - On macOS / Linux, uses `which <command>`.
 *
 * Uses `execFile` (not `exec`) so that user-provided paths containing
 * spaces or special characters are never passed through a shell.
 *
 * If `commandOrPath` looks like an absolute or relative path (contains
 * `path.sep` or starts with `.`), the function checks whether the file
 * exists using `fs.access` instead of invoking an external tool.
 *
 * @param commandOrPath - Executable name on PATH, or an absolute/relative path.
 * @returns `true` when the executable is reachable, `false` otherwise.
 */
export async function isExecutableAvailable(commandOrPath: string): Promise<boolean> {
  // If the string looks like a path (contains separator or starts with dot),
  // skip `which`/`where` and check the file directly.
  if (isFilePathLike(commandOrPath)) {
    try {
      const fs = await import('fs');
      const { promisify } = await import('util');
      const access = promisify(fs.access);
      await access(commandOrPath, fs.constants.X_OK);
      return true;
    } catch {
      return false;
    }
  }

  const cmd = process.platform === 'win32' ? 'where' : 'which';
  const { promisify } = await import('util');
  const execFileAsync = promisify(execFile);

  try {
    await execFileAsync(cmd, [commandOrPath]);
    return true;
  } catch {
    return false;
  }
}

function isFilePathLike(commandOrPath: string): boolean {
  return (
    path.isAbsolute(commandOrPath) ||
    path.win32.isAbsolute(commandOrPath) ||
    commandOrPath.includes('/') ||
    commandOrPath.includes('\\') ||
    commandOrPath.startsWith('.')
  );
}
