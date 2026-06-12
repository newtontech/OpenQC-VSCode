/**
 * Output/Log Document Detection for Standalone LSP Diagnostics
 *
 * Detects whether a document is an input file, output file, or runtime log
 * for routing diagnostics to the correct LSP server. This enables OpenQC to
 * handle output/log files from quantum chemistry codes without requiring the
 * input-language LSP to process them directly.
 *
 * @module smoke/documentDetector
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/159
 */

import type { DocumentDetectionResult, DocumentKind } from './types';
import { listBundledLspServers } from '../lsp/registry';
import type { LSPServerRegistryEntry } from '../lsp/types';

// ---------------------------------------------------------------------------
// Known output file name patterns
//
// IMPORTANT: Patterns must be specific enough to avoid false positives.
// Avoid overly broad patterns like /^.*\.out$/i which match every .out file.
// ---------------------------------------------------------------------------

/** Map of language ID to known output file name patterns. */
const OUTPUT_FILE_PATTERNS: ReadonlyMap<string, readonly RegExp[]> = new Map([
  ['vasp', [/^DOSCAR$/i, /^PROCAR$/i, /^CHG$/i, /^WAVECAR$/i, /^XDATCAR$/i, /^LOCPOT$/i]],
  ['gaussian', [/^Gaussian.*\.log$/i]],
  ['orca', [/^ORCA.*\.out$/i]],
  ['qe', [/^pw\.out$/i, /^ph\.out$/i, /^cppp\.out$/i]],
  ['gamess', [/^GAMESS.*\.log$/i, /^GAMESS.*\.dat$/i]],
  ['nwchem', [/^.*\.nwout$/i]],
  ['cp2k', [/^.*\.ener$/i, /^.*\.force$/i, /^.*-pos-1\.xyz$/i, /^.*\.restart$/i]],
  ['abacus', [/^running_.*\.log$/i, /^OUT\.ABACUS$/i]],
  ['lammps', [/^dump\..+$/i]],
  ['gromacs', [/^.*\.xvg$/i, /^.*\.trr$/i, /^.*\.xtc$/i, /^.*\.edr$/i]],
  ['gpumd', [/^thermo\.out$/i, /^xyz\.out$/i, /^velo\.out$/i]],
]);

/** Map of language ID to known log file name patterns. */
const LOG_FILE_PATTERNS: ReadonlyMap<string, readonly RegExp[]> = new Map([
  ['gaussian', [/^Gaussian.*\.log$/i]],
  ['orca', [/^ORCA.*\.out$/i]],
  ['qe', [/^pw\.out$/i]],
  ['gamess', [/^GAMESS.*\.log$/i]],
  ['nwchem', [/^.*\.nwout$/i]],
  ['cp2k', [/^cp2k.*\.log$/i]],
  ['gromacs', [/^md\.log$/i]],
  ['gpumd', [/^run\.in\.log$/i]],
]);

// ---------------------------------------------------------------------------
// Content-based heuristics
//
// These must be highly specific to avoid false matches. Use multi-word
// phrases that are unique to each code's output format.
// ---------------------------------------------------------------------------

/** Map of language ID to output file content signatures. */
const OUTPUT_CONTENT_SIGNATURES: ReadonlyMap<string, readonly RegExp[]> = new Map([
  ['vasp', [/FREE ENERGIE OF THE ION-ELECTRON SYSTEM/i, /Voluntary context switches/i]],
  ['gaussian', [/Entering Gaussian System/i, /Normal termination of Gaussian/i]],
  ['orca', [/^\* O   R   C   A \*/m, /ORCA TERMINATED NORMALLY/i]],
  ['qe', [/PROGRAM PWSCF ENDED/i, /PROGRAM PHONON ENDED/i]],
  ['gamess', [/GAMESS VERSION/i]],
  ['nwchem', [/NWChem SCF/i]],
  ['cp2k', [/^CP2K\|/m, /^ENERGY\|/m]],
  ['lammps', [/^LAMMPS\s+\(/m]],
  ['gromacs', [/^GROMACS\s+/m]],
]);

// ---------------------------------------------------------------------------
// Detection engine
// ---------------------------------------------------------------------------

/**
 * Detect the document kind for a given file path and optional content.
 *
 * Uses file name patterns first, then falls back to content-based detection
 * when content is provided.
 *
 * @param fileName - The base file name (not the full path).
 * @param content - Optional file content for content-based detection.
 * @returns Detection result with kind, associated server, and confidence.
 */
export function detectDocument(fileName: string, content?: string): DocumentDetectionResult {
  const servers = listBundledLspServers();

  // Phase 1: Check against input file patterns from the registry.
  // Uses two-pass matching: exact fileName matches first (higher priority),
  // then extension matches. This ensures 'run.in' matches GPUMD's exact
  // fileName before QE's '.in' extension.
  const inputMatch = matchInputFile(fileName, servers);
  if (inputMatch) {
    return {
      kind: 'input',
      serverId: inputMatch.id,
      languageId: inputMatch.languageId,
      confidence: 1.0,
    };
  }

  // Phase 2: Check output file name patterns
  const outputNameMatch = matchOutputFileName(fileName, servers);
  if (outputNameMatch) {
    return {
      kind: 'output',
      serverId: outputNameMatch.serverId,
      languageId: outputNameMatch.languageId,
      confidence: 0.9,
    };
  }

  // Phase 3: Check log file name patterns
  const logNameMatch = matchLogFileName(fileName, servers);
  if (logNameMatch) {
    return {
      kind: 'log',
      serverId: logNameMatch.serverId,
      languageId: logNameMatch.languageId,
      confidence: 0.85,
    };
  }

  // Phase 4: Content-based detection (if content provided)
  if (content !== undefined) {
    const contentMatch = matchContent(content, servers);
    if (contentMatch) {
      return {
        kind: contentMatch.kind,
        serverId: contentMatch.serverId,
        languageId: contentMatch.languageId,
        confidence: contentMatch.confidence,
      };
    }
  }

  return {
    kind: 'unknown',
    confidence: 0,
  };
}

// ---------------------------------------------------------------------------
// Internal matchers
// ---------------------------------------------------------------------------

/**
 * Match input files using a two-pass approach:
 * 1. First pass: check exact fileName matches across ALL servers.
 * 2. Second pass: check extension matches across all servers.
 *
 * This ensures that exact fileName matches (e.g., GPUMD's "run.in") take
 * priority over extension matches (e.g., QE's ".in").
 */
function matchInputFile(
  fileName: string,
  servers: readonly LSPServerRegistryEntry[]
): LSPServerRegistryEntry | undefined {
  const baseName = fileName.toLowerCase();

  // Pass 1: Exact fileName matches (highest priority)
  for (const server of servers) {
    if (server.fileNames.some(fn => fn.toLowerCase() === baseName)) {
      return server;
    }
  }

  // Pass 2: Extension matches
  for (const server of servers) {
    for (const ext of server.fileExtensions) {
      if (baseName.endsWith('.' + ext.toLowerCase())) {
        return server;
      }
    }
  }

  return undefined;
}

function matchOutputFileName(
  fileName: string,
  _servers: readonly LSPServerRegistryEntry[]
): { readonly serverId: string; readonly languageId: string } | undefined {
  for (const [languageId, patterns] of OUTPUT_FILE_PATTERNS) {
    for (const pattern of patterns) {
      if (pattern.test(fileName)) {
        const server = _servers.find(s => s.languageId === languageId);
        if (server) {
          return { serverId: server.id, languageId };
        }
      }
    }
  }

  return undefined;
}

function matchLogFileName(
  fileName: string,
  _servers: readonly LSPServerRegistryEntry[]
): { readonly serverId: string; readonly languageId: string } | undefined {
  for (const [languageId, patterns] of LOG_FILE_PATTERNS) {
    for (const pattern of patterns) {
      if (pattern.test(fileName)) {
        const server = _servers.find(s => s.languageId === languageId);
        if (server) {
          return { serverId: server.id, languageId };
        }
      }
    }
  }

  return undefined;
}

function matchContent(
  content: string,
  _servers: readonly LSPServerRegistryEntry[]
):
  | {
      readonly kind: DocumentKind;
      readonly serverId: string;
      readonly languageId: string;
      readonly confidence: number;
    }
  | undefined {
  for (const [languageId, signatures] of OUTPUT_CONTENT_SIGNATURES) {
    for (const signature of signatures) {
      if (signature.test(content)) {
        const server = _servers.find(s => s.languageId === languageId);
        if (server) {
          return {
            kind: 'output',
            serverId: server.id,
            languageId,
            confidence: 0.8,
          };
        }
      }
    }
  }

  return undefined;
}

// ---------------------------------------------------------------------------
// Utility
// ---------------------------------------------------------------------------

/**
 * Check if a file name looks like a known output or log file (without content).
 *
 * @param fileName - The base file name.
 * @returns True if the file name matches known output/log patterns.
 */
export function isOutputOrLogFile(fileName: string): boolean {
  const result = detectDocument(fileName);
  return result.kind === 'output' || result.kind === 'log';
}

/**
 * Get all known output file patterns for a given language ID.
 *
 * @param languageId - The VS Code language ID.
 * @returns Array of RegExp patterns for output files of that language.
 */
export function getOutputPatternsForLanguage(languageId: string): readonly RegExp[] {
  return OUTPUT_FILE_PATTERNS.get(languageId) ?? [];
}

/**
 * Get all known log file patterns for a given language ID.
 *
 * @param languageId - The VS Code language ID.
 * @returns Array of RegExp patterns for log files of that language.
 */
export function getLogPatternsForLanguage(languageId: string): readonly RegExp[] {
  return LOG_FILE_PATTERNS.get(languageId) ?? [];
}
