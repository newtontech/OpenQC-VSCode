import * as path from 'path';
import {
  getCapabilityManifestPath,
  loadCapabilityManifest,
  validateCapabilityManifest,
} from '../../../src/lsp/capabilityManifest';
import { getLspServerByLanguageId } from '../../../src/lsp/registry';

const REPO_ROOT = path.resolve(__dirname, '../../../');

describe('capabilityManifest', () => {
  it('loads and validates a sibling LSP capability manifest', () => {
    const server = getLspServerByLanguageId('dpgen');
    expect(server).toBeDefined();

    const result = loadCapabilityManifest(REPO_ROOT, server!);
    expect(result.error).toBeUndefined();
    expect(result.manifest).toBeDefined();
    expect(result.manifest!.id).toBe('dpgen-lsp');
    expect(result.manifest!.agentCli.operations).toEqual(
      expect.arrayContaining([
        'capabilities',
        'check',
        'context',
        'complete',
        'hover',
        'symbols',
        'fix',
      ])
    );
  });

  it('returns the expected sibling manifest path', () => {
    const server = getLspServerByLanguageId('vasp');
    expect(server).toBeDefined();
    expect(getCapabilityManifestPath(REPO_ROOT, server!)).toContain(
      'VASP-LSP/lsp-capabilities.json'
    );
  });

  it('reports registry mismatches', () => {
    const server = getLspServerByLanguageId('qe');
    expect(server).toBeDefined();
    const error = validateCapabilityManifest(
      {
        schema: 'OpenQCLspCapabilities',
        version: 1,
        id: 'wrong',
        languageId: 'qe',
        executable: 'qe-lsp',
        defaultBranch: 'main',
        filePatterns: ['*.in'],
        blockingPolicy: { mode: 'warning-only' },
        capabilities: ['diagnostics'],
        agentCli: {
          command: 'qe-lsp-tool',
          operations: ['capabilities', 'check', 'context', 'complete', 'hover', 'symbols', 'fix'],
          jsonFormat: true,
          failOnBlocking: true,
        },
        diagnosticSchema: 'diagnostics/diagnostic-engine-v1.schema.json',
        wikiPaths: {},
        fixturePaths: {},
        outputLogPatterns: [],
        openqc: {
          registryId: 'qe-lsp',
          repoName: 'qe-lsp',
          contextContract: 'DSLAuthoringContext',
          diagnosticEnvelope: 'DiagnosticEnvelope/v1',
        },
      },
      server!
    );
    expect(error).toContain('id wrong != qe-lsp');
  });
});
