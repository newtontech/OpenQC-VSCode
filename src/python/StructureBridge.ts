/**
 * TypeScript wrapper for the Python structure bridge.
 * @module python/StructureBridge
 */

import { execPythonJson, type BridgeResponse } from './PythonBridge';
import type { OpenQCStructure } from '../structures/OpenQCStructure';

export async function parseStructure(
  filePath: string,
  formatHint?: string,
  content?: string
): Promise<BridgeResponse<OpenQCStructure>> {
  return execPythonJson<OpenQCStructure>('openqc.bridge.structure_bridge', {
    command: 'parse',
    args: { path: filePath, format: formatHint || 'auto', content },
  });
}

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
