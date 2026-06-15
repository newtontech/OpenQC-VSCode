import * as fs from 'fs';
import * as path from 'path';
import type { LspCapabilityManifest, LSPServerRegistryEntry } from './types';

export interface CapabilityManifestLoadResult {
  readonly manifest?: LspCapabilityManifest;
  readonly manifestPath: string;
  readonly error?: string;
}

const REQUIRED_AGENT_OPERATIONS = [
  'capabilities',
  'check',
  'context',
  'complete',
  'hover',
  'symbols',
  'fix',
] as const;

const REQUIRED_CAPABILITIES = [
  'diagnostics',
  'rich-diagnostics',
  'completion',
  'hover',
  'symbols',
  'fix-preview',
  'llm-wiki',
  'openqc-context',
] as const;

const REQUIRED_WIKI_PATHS = [
  'plan',
  'rawAssets',
  'entities',
  'concepts',
  'synthesis',
  'index',
  'log',
] as const;

const REQUIRED_FIXTURE_PATHS = ['valid', 'invalid', 'logs'] as const;

export function getSiblingRepoPath(openqcRepoRoot: string, server: LSPServerRegistryEntry): string {
  const repoName = server.localLaunch?.repoName ?? server.repository.split('/').pop() ?? server.id;
  return path.resolve(openqcRepoRoot, '..', repoName);
}

export function getCapabilityManifestPath(
  openqcRepoRoot: string,
  server: LSPServerRegistryEntry
): string {
  return path.join(getSiblingRepoPath(openqcRepoRoot, server), 'lsp-capabilities.json');
}

export function loadCapabilityManifest(
  openqcRepoRoot: string,
  server: LSPServerRegistryEntry
): CapabilityManifestLoadResult {
  const manifestPath = getCapabilityManifestPath(openqcRepoRoot, server);
  try {
    const raw = fs.readFileSync(manifestPath, 'utf8');
    const parsed = JSON.parse(raw) as unknown;
    const validationError = validateCapabilityManifest(parsed, server);
    if (validationError) {
      return { manifestPath, error: validationError };
    }
    return { manifest: parsed as LspCapabilityManifest, manifestPath };
  } catch (error) {
    return {
      manifestPath,
      error: error instanceof Error ? error.message : String(error),
    };
  }
}

export function validateCapabilityManifest(
  raw: unknown,
  server: LSPServerRegistryEntry
): string | undefined {
  if (!raw || typeof raw !== 'object') {
    return 'manifest is not an object';
  }
  const manifest = raw as Partial<LspCapabilityManifest>;
  const mismatches: string[] = [];

  if (manifest.schema !== 'OpenQCLspCapabilities') {
    mismatches.push('schema must be OpenQCLspCapabilities');
  }
  if (manifest.version !== 1) {
    mismatches.push('version must be 1');
  }
  if (manifest.id !== server.id) {
    mismatches.push(`id ${String(manifest.id)} != ${server.id}`);
  }
  if (manifest.languageId !== server.languageId) {
    mismatches.push(`languageId ${String(manifest.languageId)} != ${server.languageId}`);
  }
  if (manifest.executable !== server.executable) {
    mismatches.push(`executable ${String(manifest.executable)} != ${server.executable}`);
  }
  if (manifest.defaultBranch !== server.defaultBranch) {
    mismatches.push(`defaultBranch ${String(manifest.defaultBranch)} != ${server.defaultBranch}`);
  }
  if (!Array.isArray(manifest.filePatterns) || manifest.filePatterns.length === 0) {
    mismatches.push('filePatterns must be a non-empty array');
  }
  if (!manifest.blockingPolicy || typeof manifest.blockingPolicy !== 'object') {
    mismatches.push('blockingPolicy is required');
  } else if (!['blocking', 'warning-only'].includes(String(manifest.blockingPolicy.mode))) {
    mismatches.push('blockingPolicy.mode must be blocking or warning-only');
  }
  if (!Array.isArray(manifest.capabilities) || manifest.capabilities.length === 0) {
    mismatches.push('capabilities must be a non-empty array');
  } else {
    const capabilities = new Set(manifest.capabilities);
    for (const capability of REQUIRED_CAPABILITIES) {
      if (!capabilities.has(capability)) {
        mismatches.push(`capabilities missing ${capability}`);
      }
    }
  }
  if (!manifest.agentCli || typeof manifest.agentCli !== 'object') {
    mismatches.push('agentCli is required');
  } else {
    if (!manifest.agentCli.command) {
      mismatches.push('agentCli.command is required');
    }
    const operations = new Set(manifest.agentCli.operations ?? []);
    for (const operation of REQUIRED_AGENT_OPERATIONS) {
      if (!operations.has(operation)) {
        mismatches.push(`agentCli.operations missing ${operation}`);
      }
    }
    if (manifest.agentCli.jsonFormat !== true) {
      mismatches.push('agentCli.jsonFormat must be true');
    }
    if (manifest.agentCli.failOnBlocking !== true) {
      mismatches.push('agentCli.failOnBlocking must be true');
    }
  }
  if (!manifest.diagnosticSchema) {
    mismatches.push('diagnosticSchema is required');
  }
  if (!manifest.wikiPaths) {
    mismatches.push('wikiPaths is required');
  } else {
    for (const wikiPath of REQUIRED_WIKI_PATHS) {
      if (!manifest.wikiPaths[wikiPath]) {
        mismatches.push(`wikiPaths.${wikiPath} is required`);
      }
    }
  }
  if (!manifest.fixturePaths) {
    mismatches.push('fixturePaths is required');
  } else {
    for (const fixturePath of REQUIRED_FIXTURE_PATHS) {
      if (!Array.isArray(manifest.fixturePaths[fixturePath])) {
        mismatches.push(`fixturePaths.${fixturePath} must be an array`);
      }
    }
  }
  if (!manifest.openqc || manifest.openqc.registryId !== server.id) {
    mismatches.push('openqc.registryId must match registry id');
  } else {
    if (manifest.openqc.contextContract !== 'DSLAuthoringContext') {
      mismatches.push('openqc.contextContract must be DSLAuthoringContext');
    }
    if (manifest.openqc.diagnosticEnvelope !== 'DiagnosticEnvelope/v1') {
      mismatches.push('openqc.diagnosticEnvelope must be DiagnosticEnvelope/v1');
    }
  }

  return mismatches.length ? mismatches.join('; ') : undefined;
}
