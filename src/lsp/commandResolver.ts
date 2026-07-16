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
import * as fs from 'fs';
import * as os from 'os';
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
  | {
      kind: 'pathOrCommand';
      command: string;
      args: string[];
      env?: Record<string, string>;
      cwd?: string;
    }
  | {
      kind: 'pythonModule';
      python: string;
      module: string;
      args: string[];
      env?: Record<string, string>;
      cwd?: string;
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

export interface LspCommandResolutionOptions {
  /** Extension install/source path used to locate sibling LSP repositories. */
  extensionPath?: string;
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
  overrides: LspCommandOverrides,
  options: LspCommandResolutionOptions = {}
): ResolvedLspCommand {
  if (!hasUserOverride(overrides)) {
    const localCommand = resolveLocalRepositoryCommand(entry, overrides, options);
    if (localCommand) {
      return localCommand;
    }
  }

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
 *
 * Only user/workspace/language overrides are returned. Extension-contributed
 * defaults from package.json are intentionally ignored so local sibling LSP
 * checkouts can still be auto-detected.
 */
export function readCommandOverrides(
  config: vscode.WorkspaceConfiguration,
  languageId: string
): LspCommandOverrides {
  return {
    path: readExplicitConfigValue<string>(config, `${languageId}.path`),
    command: readExplicitConfigValue<string>(config, `${languageId}.command`),
    args: readExplicitConfigValue<string[]>(config, `${languageId}.args`),
    env: readExplicitConfigValue<Record<string, string>>(config, `${languageId}.env`),
  };
}

function readExplicitConfigValue<T>(
  config: vscode.WorkspaceConfiguration,
  key: string
): T | undefined {
  const inspection = typeof config.inspect === 'function' ? config.inspect<T>(key) : undefined;
  if (inspection) {
    return firstDefined<T>(
      inspection.workspaceFolderLanguageValue,
      inspection.workspaceLanguageValue,
      inspection.globalLanguageValue,
      inspection.workspaceFolderValue,
      inspection.workspaceValue,
      inspection.globalValue
    );
  }

  return config.get<T | undefined>(key, undefined);
}

function firstDefined<T>(...values: Array<T | undefined>): T | undefined {
  return values.find((value): value is T => value !== undefined);
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

function hasUserOverride(overrides: LspCommandOverrides): boolean {
  return (
    overrides.command !== undefined ||
    overrides.path !== undefined ||
    overrides.args !== undefined ||
    overrides.env !== undefined
  );
}

function resolveLocalRepositoryCommand(
  entry: LSPServerRegistryEntry,
  overrides: LspCommandOverrides,
  options: LspCommandResolutionOptions
): ResolvedLspCommand | undefined {
  if (!entry.localLaunch) {
    return undefined;
  }

  const repoRoot = findSiblingRepository(entry.localLaunch.repoName, options.extensionPath);
  if (!repoRoot) {
    return undefined;
  }

  const args = overrides.args ?? [...(entry.args ?? ['--stdio'])];
  const env = overrides.env || {};

  if (entry.localLaunch.kind === 'pythonFunction') {
    const sourcePath = path.join(repoRoot, entry.localLaunch.sourcePath || 'src');
    const pythonPath = mergePathList(sourcePath, process.env.PYTHONPATH);
    const code = [
      `from ${entry.localLaunch.importPath} import ${entry.localLaunch.functionName} as _main`,
      'raise SystemExit(_main())',
    ].join('; ');

    return {
      kind: 'pathOrCommand',
      command: process.env.OPENQC_LSP_PYTHON || 'python3',
      args: ['-c', code, ...args],
      cwd: repoRoot,
      env: { ...env, PYTHONPATH: pythonPath },
    };
  }

  if (entry.localLaunch.kind === 'nodeScript') {
    const script = path.join(repoRoot, entry.localLaunch.scriptPath);
    if (!fs.existsSync(script)) {
      return undefined;
    }
    return {
      kind: 'pathOrCommand',
      command: process.env.OPENQC_LSP_NODE || 'node',
      args: [script, ...args],
      cwd: repoRoot,
      env: Object.keys(env).length > 0 ? env : undefined,
    };
  }

  const binary = path.join(repoRoot, 'target', 'debug', entry.localLaunch.binaryName);
  if (fs.existsSync(binary)) {
    return {
      kind: 'pathOrCommand',
      command: binary,
      args,
      cwd: repoRoot,
      env: Object.keys(env).length > 0 ? env : undefined,
    };
  }

  return {
    kind: 'pathOrCommand',
    command: process.env.OPENQC_LSP_CARGO || 'cargo',
    args: [
      'run',
      '--quiet',
      '--bin',
      entry.localLaunch.cargoBin || entry.localLaunch.binaryName,
      '--',
    ],
    cwd: repoRoot,
    env: Object.keys(env).length > 0 ? env : undefined,
  };
}

function findSiblingRepository(repoName: string, extensionPath?: string): string | undefined {
  const roots = candidateSearchRoots(extensionPath);
  for (const root of roots) {
    const candidate = path.join(root, repoName);
    if (fs.existsSync(candidate) && fs.statSync(candidate).isDirectory()) {
      return candidate;
    }
  }
  return undefined;
}

function candidateSearchRoots(extensionPath?: string): string[] {
  const roots: string[] = [];
  const add = (candidate?: string) => {
    if (!candidate) {
      return;
    }
    const normalized = path.resolve(candidate);
    if (!roots.includes(normalized)) {
      roots.push(normalized);
    }
  };

  add(process.env.OPENQC_LSP_REPOSITORY_ROOT);

  let current = path.resolve(extensionPath || process.cwd());
  for (let depth = 0; depth < 5; depth += 1) {
    add(path.join(current, '.worktrees-lsp-latest'));
    add(path.join(path.dirname(current), '.worktrees-lsp-latest'));
    add(current);
    add(path.dirname(current));
    current = path.dirname(current);
  }

  add(path.join(os.homedir(), 'Desktop', 'code', '.worktrees-lsp-latest'));
  add(path.join(os.homedir(), 'Desktop', 'code'));
  return roots;
}

function mergePathList(first: string, existing?: string): string {
  return existing && existing.length > 0 ? `${first}${path.delimiter}${existing}` : first;
}
