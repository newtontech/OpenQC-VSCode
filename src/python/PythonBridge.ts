/**
 * PythonBridge - Robust subprocess bridge for spawning Python commands.
 *
 * Provides:
 * - JSON stdin/stdout protocol.
 * - Request timeout and cancellation.
 * - Deterministic error handling for missing executable, invalid JSON,
 *   subprocess non-zero exit, and timeout cases.
 *
 * @module python/PythonBridge
 */

import { execFile } from 'child_process';
import * as path from 'path';
import * as vscode from 'vscode';
import type { BackendCheckResult, BridgeError } from './PythonBackendStatus';
import { Logger } from '../utils/Logger';

const logger = Logger.getInstance();

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

const DEFAULT_TIMEOUT_MS = 30_000;

function getPythonPath(): string {
  const config = vscode.workspace.getConfiguration('openqc');
  return config.get<string>('pythonPath', 'python3');
}

function getTimeoutMs(): number {
  const config = vscode.workspace.getConfiguration('openqc');
  return config.get<number>('python.timeoutMs', DEFAULT_TIMEOUT_MS);
}

function getBundledPythonPath(): string {
  return path.resolve(__dirname, '..', '..', 'python');
}

function getPythonEnv(): NodeJS.ProcessEnv {
  const bundledPythonPath = getBundledPythonPath();
  const existing = process.env.PYTHONPATH;
  return {
    ...process.env,
    PYTHONPATH: existing ? `${bundledPythonPath}${path.delimiter}${existing}` : bundledPythonPath,
  };
}

// ---------------------------------------------------------------------------
// Bridge response
// ---------------------------------------------------------------------------

export interface BridgeResponse<T> {
  success: boolean;
  data?: T;
  error?: BridgeError;
  timedOut?: boolean;
  exitCode?: number | null;
  stderr?: string;
}

// ---------------------------------------------------------------------------
// Core exec helper
// ---------------------------------------------------------------------------

/**
 * Execute a Python module with a JSON payload on stdin and return parsed
 * JSON from stdout.
 *
 * @param modulePath - Python module path (e.g. "openqc.bridge.check_backend").
 * @param payload    - JSON-serializable object to send on stdin.
 * @param options    - Optional overrides for timeout and cancellation.
 */
export function execPythonJson<T = unknown>(
  modulePath: string,
  payload?: Record<string, unknown>,
  options?: {
    timeoutMs?: number;
    cancelToken?: vscode.CancellationToken;
    pythonPath?: string;
  }
): Promise<BridgeResponse<T>> {
  const pythonPath = options?.pythonPath ?? getPythonPath();
  const timeoutMs = options?.timeoutMs ?? getTimeoutMs();

  return new Promise<BridgeResponse<T>>(resolve => {
    const args = ['-m', modulePath];
    const stdinData = payload ? JSON.stringify(payload) : undefined;

    logger.debug(`PythonBridge: spawning ${pythonPath} ${args.join(' ')}`);

    const child = execFile(
      pythonPath,
      args,
      {
        timeout: timeoutMs,
        maxBuffer: 10 * 1024 * 1024, // 10 MB
        env: getPythonEnv(),
      },
      (error, stdout, stderr) => {
        if (options?.cancelToken?.isCancellationRequested) {
          resolve({
            success: false,
            error: { message: 'Request was cancelled' },
            timedOut: false,
          });
          return;
        }

        if (error) {
          // Timeout case
          if ((error as any).killed) {
            resolve({
              success: false,
              error: { message: `Python process timed out after ${timeoutMs}ms` },
              timedOut: true,
              exitCode: null,
              stderr: stderr?.trim() || undefined,
            });
            return;
          }

          // Missing executable
          if ((error as any).code === 'ENOENT') {
            resolve({
              success: false,
              error: {
                message: `Python executable not found: ${pythonPath}`,
                hint: 'Install Python 3.8+ or set openqc.pythonPath in settings.',
              },
              exitCode: null,
              stderr: stderr?.trim() || undefined,
            });
            return;
          }

          // Non-zero exit
          resolve({
            success: false,
            error: { message: error.message },
            exitCode: typeof error.code === 'number' ? error.code : null,
            stderr: stderr?.trim() || undefined,
          });
          return;
        }

        // Try to parse stdout as JSON
        const trimmed = stdout.trim();
        if (!trimmed) {
          resolve({
            success: false,
            error: { message: 'Python process returned empty output' },
            stderr: stderr?.trim() || undefined,
          });
          return;
        }

        try {
          const data = JSON.parse(trimmed) as T;
          resolve({ success: true, data });
        } catch (parseError) {
          // Check if stderr has structured error
          const stderrTrimmed = stderr?.trim() || '';
          let bridgeError: BridgeError = {
            message: `Invalid JSON response from Python: ${parseError}`,
            rawStdout: trimmed.substring(0, 500),
          };

          if (stderrTrimmed) {
            try {
              const stderrJson = JSON.parse(stderrTrimmed);
              if (stderrJson.error) {
                bridgeError = stderrJson.error;
              }
            } catch {
              bridgeError.rawStderr = stderrTrimmed.substring(0, 500);
            }
          }

          resolve({
            success: false,
            error: bridgeError,
            stderr: stderrTrimmed || undefined,
          });
        }
      }
    );

    // Send stdin data
    if (stdinData && child.stdin) {
      child.stdin.write(stdinData);
      child.stdin.end();
    }

    // Handle cancellation
    if (options?.cancelToken) {
      options.cancelToken.onCancellationRequested(() => {
        child.kill();
      });
    }
  });
}

// ---------------------------------------------------------------------------
// Backend health check
// ---------------------------------------------------------------------------

/**
 * Check the scientific Python backend status.
 *
 * Returns Python version, available packages, and external tools.
 */
export async function checkBackend(options?: {
  timeoutMs?: number;
  cancelToken?: vscode.CancellationToken;
  pythonPath?: string;
}): Promise<BridgeResponse<BackendCheckResult>> {
  return execPythonJson<BackendCheckResult>('openqc.bridge.check_backend', undefined, options);
}
