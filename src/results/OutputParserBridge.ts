/**
 * TypeScript wrapper for the cclib-powered output parser bridge.
 *
 * @module results/OutputParserBridge
 */

import { execPythonJson, type BridgeResponse } from '../python/PythonBridge';
import type { OpenQCResults } from './OpenQCResults';
import { Logger } from '../utils/Logger';

const logger = Logger.getInstance();

export interface OutputTrajectory {
  sourceFile?: string;
  software?: string;
  supported: boolean;
  frames: unknown[];
  frameCount: number;
  energies?: number[];
  cclibAvailable?: boolean;
  warnings?: string[];
}

/**
 * Parse a calculation output file.
 *
 * Uses cclib when available, with native fallback for basic extraction.
 *
 * @param filePath - Path to the output file.
 * @param software - Optional software hint (e.g. 'gaussian', 'orca').
 */
export async function parseOutput(
  filePath: string,
  software?: string
): Promise<BridgeResponse<OpenQCResults>> {
  const args: Record<string, unknown> = {
    path: filePath,
    software: software || 'auto',
  };

  const result = await execPythonJson<OpenQCResults>('openqc.bridge.output_bridge', {
    command: 'parse',
    args,
  });

  if (result.success && result.data) {
    logger.info(`Output parsed: ${filePath}, software=${result.data.software}`);

    if (result.data.cclibAvailable === false) {
      logger.warn('cclib not available — output parsing used native fallback');
    }
  }

  return result;
}

/**
 * Get a brief summary of calculation results.
 */
export async function summarizeOutput(
  filePath: string,
  software?: string
): Promise<BridgeResponse<Record<string, unknown>>> {
  return execPythonJson<Record<string, unknown>>('openqc.bridge.output_bridge', {
    command: 'summarize',
    args: { path: filePath, software: software || 'auto' },
  });
}

/**
 * Extract an optimization trajectory from a calculation output file.
 *
 * Returns a structured unsupported response when the output format or parser
 * does not expose coordinate frames.
 */
export async function extractTrajectory(
  filePath: string,
  software?: string
): Promise<BridgeResponse<OutputTrajectory>> {
  return execPythonJson<OutputTrajectory>('openqc.bridge.output_bridge', {
    command: 'trajectory',
    args: { path: filePath, software: software || 'auto' },
  });
}
