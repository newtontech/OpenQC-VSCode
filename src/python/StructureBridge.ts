/**
 * TypeScript wrapper for the Python structure bridge.
 *
 * Uses the shared PythonBridge subprocess module to invoke
 * structure parsing, conversion, and supercell generation.
 *
 * @module python/StructureBridge
 */

import { execPythonJson, type BridgeResponse } from './PythonBridge';
import type { OpenQCStructure } from '../structures/OpenQCStructure';
import { validateOpenQCStructure } from '../structures/validation';
import { Logger } from '../utils/Logger';

const logger = Logger.getInstance();

// ---------------------------------------------------------------------------
// Commands
// ---------------------------------------------------------------------------

/**
 * Parse a structure file using the Python bridge.
 *
 * Tries native parsers first, then falls back to pymatgen/ASE.
 *
 * @param filePath - Path to the structure file.
 * @param formatHint - Optional format override (e.g. 'cif', 'poscar', 'xyz').
 * @param content - Optional file content (avoids reading from disk in native parsers).
 */
export async function parseStructure(
  filePath: string,
  formatHint?: string,
  content?: string
): Promise<BridgeResponse<OpenQCStructure>> {
  const args: Record<string, unknown> = {
    path: filePath,
    format: formatHint || 'auto',
  };
  if (content) {
    args.content = content;
  }

  const result = await execPythonJson<OpenQCStructure>('openqc.bridge.structure_bridge', {
    command: 'parse',
    args,
  });

  if (result.success && result.data) {
    const validation = validateOpenQCStructure(result.data);
    if (!validation.valid) {
      logger.warn(`Parsed structure validation warnings: ${validation.errors.join('; ')}`);
    }
  }

  return result;
}

/**
 * Generate a supercell from an existing OpenQCStructure.
 *
 * @param structure - The base structure (must have cell).
 * @param na - Supercell multiplier along a.
 * @param nb - Supercell multiplier along b.
 * @param nc - Supercell multiplier along c.
 */
export async function generateSupercell(
  structure: OpenQCStructure,
  na: number = 2,
  nb: number = 2,
  nc: number = 2
): Promise<BridgeResponse<OpenQCStructure>> {
  return execPythonJson<OpenQCStructure>('openqc.bridge.structure_bridge', {
    command: 'supercell',
    args: { structure, na, nb, nc },
  });
}

/**
 * Convert a structure file to another format via the Python bridge.
 */
export async function convertStructure(
  filePath: string,
  targetFormat: string,
  formatHint?: string
): Promise<BridgeResponse<OpenQCStructure>> {
  return execPythonJson<OpenQCStructure>('openqc.bridge.structure_bridge', {
    command: 'convert',
    args: { path: filePath, to: targetFormat, format: formatHint || 'auto' },
  });
}
