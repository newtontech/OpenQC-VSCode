/**
 * OpenQC Tool Registry - Exposes operations as VS Code Language Model Tools
 * and an optional MCP server.
 *
 * Tool classification:
 * - Read-only: No confirmation needed (parse, explain, check)
 * - Preview-only: No source modification (preview, supercell, plot)
 * - Write: Requires diff/output confirmation (convert, generate, export)
 * - External: Requires explicit confirmation (Multiwfn, c2x)
 *
 * @module agent/OpenQCToolRegistry
 */

// ---------------------------------------------------------------------------
// Tool definitions
// ---------------------------------------------------------------------------

export type ToolCategory = 'read-only' | 'preview-only' | 'write' | 'external';

export interface OpenQCToolDefinition {
  name: string;
  description: string;
  category: ToolCategory;
  inputSchema: Record<string, unknown>;
  handler: string; // function reference name
}

export const OPENQC_TOOLS: OpenQCToolDefinition[] = [
  // Read-only tools
  {
    name: 'openqc.parseStructure',
    description: 'Parse a molecular/crystal structure file and return OpenQCStructure JSON',
    category: 'read-only',
    inputSchema: {
      type: 'object',
      properties: {
        path: { type: 'string', description: 'File path' },
        format: {
          type: 'string',
          description: 'Format hint (auto, xyz, poscar, cif, gaussian, orca)',
        },
      },
      required: ['path'],
    },
    handler: 'parseStructure',
  },
  {
    name: 'openqc.parseCalculationOutput',
    description: 'Parse a quantum chemistry calculation output file using cclib',
    category: 'read-only',
    inputSchema: {
      type: 'object',
      properties: {
        path: { type: 'string', description: 'Output file path' },
        software: { type: 'string', description: 'Software hint (auto, gaussian, orca, nwchem)' },
      },
      required: ['path'],
    },
    handler: 'parseOutput',
  },
  {
    name: 'openqc.checkPythonBackend',
    description: 'Check Scientific Python backend: Python version, packages, external tools',
    category: 'read-only',
    inputSchema: { type: 'object', properties: {} },
    handler: 'checkBackend',
  },
  {
    name: 'openqc.checkExternalAnalyzers',
    description: 'Check availability and configuration of external analyzers (Multiwfn, c2x)',
    category: 'read-only',
    inputSchema: { type: 'object', properties: {} },
    handler: 'checkAnalyzers',
  },
  {
    name: 'openqc.explainInput',
    description: 'Explain the parameters and keywords in a quantum chemistry input file',
    category: 'read-only',
    inputSchema: {
      type: 'object',
      properties: {
        path: { type: 'string', description: 'Input file path' },
      },
      required: ['path'],
    },
    handler: 'explainInput',
  },

  // Preview-only tools
  {
    name: 'openqc.previewStructure',
    description: 'Preview a structure in the 3D viewer without modifying the source file',
    category: 'preview-only',
    inputSchema: {
      type: 'object',
      properties: {
        path: { type: 'string', description: 'File path' },
      },
      required: ['path'],
    },
    handler: 'previewStructure',
  },
  {
    name: 'openqc.generateSupercell',
    description: 'Generate a supercell preview without modifying the source file',
    category: 'preview-only',
    inputSchema: {
      type: 'object',
      properties: {
        path: { type: 'string', description: 'Structure file path' },
        na: { type: 'number', description: 'Repeat along a' },
        nb: { type: 'number', description: 'Repeat along b' },
        nc: { type: 'number', description: 'Repeat along c' },
      },
      required: ['path'],
    },
    handler: 'generateSupercell',
  },
  {
    name: 'openqc.plotSCFConvergence',
    description: 'Plot SCF convergence from a calculation output file',
    category: 'preview-only',
    inputSchema: {
      type: 'object',
      properties: {
        path: { type: 'string', description: 'Output file path' },
      },
      required: ['path'],
    },
    handler: 'plotSCF',
  },

  // Write tools
  {
    name: 'openqc.convertStructure',
    description: 'Convert a structure to another format (requires output path confirmation)',
    category: 'write',
    inputSchema: {
      type: 'object',
      properties: {
        path: { type: 'string', description: 'Source file path' },
        to: { type: 'string', description: 'Target format (xyz, cif, poscar, gaussian)' },
        outputPath: { type: 'string', description: 'Output file path' },
      },
      required: ['path', 'to'],
    },
    handler: 'convertStructure',
  },
  {
    name: 'openqc.generateInput',
    description: 'Generate a quantum chemistry input file from a template',
    category: 'write',
    inputSchema: {
      type: 'object',
      properties: {
        software: { type: 'string', description: 'Target software (gaussian, orca, cp2k, qe)' },
        task: { type: 'string', description: 'Calculation task (optimization, frequency, sp)' },
        structurePath: { type: 'string', description: 'Structure file path' },
        method: { type: 'string', description: 'Method (B3LYP, PBE, HF)' },
        basis: { type: 'string', description: 'Basis set (6-31G*, def2-SVP)' },
      },
      required: ['software', 'task', 'structurePath'],
    },
    handler: 'generateInput',
  },
];

// ---------------------------------------------------------------------------
// Tool lookup
// ---------------------------------------------------------------------------

export function getToolByName(name: string): OpenQCToolDefinition | undefined {
  return OPENQC_TOOLS.find(t => t.name === name);
}

export function getToolsByCategory(category: ToolCategory): OpenQCToolDefinition[] {
  return OPENQC_TOOLS.filter(t => t.category === category);
}

export function getToolList(): Array<{
  name: string;
  description: string;
  category: ToolCategory;
}> {
  return OPENQC_TOOLS.map(t => ({
    name: t.name,
    description: t.description,
    category: t.category,
  }));
}

// ---------------------------------------------------------------------------
// Safety classification
// ---------------------------------------------------------------------------

export function requiresConfirmation(toolName: string): boolean {
  const tool = getToolByName(toolName);
  if (!tool) return true; // Unknown tools require confirmation
  return tool.category === 'write' || tool.category === 'external';
}

export function isReadOnly(toolName: string): boolean {
  const tool = getToolByName(toolName);
  return tool?.category === 'read-only';
}
