/**
 * OpenQC ↔ Bohrium registry alignment unit tests.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/194
 */

import { execFileSync } from 'child_process';
import * as fs from 'fs';
import * as path from 'path';
import { getLspDiagnosticReadiness, listBundledLspServers } from '../../../src/lsp/registry';

const REPO_ROOT = path.resolve(__dirname, '../../../');
const CODE_ROOT = path.resolve(REPO_ROOT, '..');
const DEFAULT_BOHRIUM_REGISTRY = path.join(
  CODE_ROOT,
  'bohrium_skills/bohrium_skills/lsp-skills/references/lsp_backends.yaml'
);
const ALIGNMENT_SCRIPT = path.join(REPO_ROOT, 'scripts/check-bohrium-registry-alignment.mjs');

function loadBohriumBackends(registryPath: string): Array<{
  id: string;
  agent_cli?: string;
}> {
  const registry = JSON.parse(fs.readFileSync(registryPath, 'utf8')) as {
    backends: Array<{ id: string; agent_cli?: string }>;
  };
  return registry.backends;
}

describe('OpenQC ↔ Bohrium registry alignment', () => {
  const openqcServers = listBundledLspServers();
  const registryPath = process.env.BOHRIUM_LSP_REGISTRY ?? DEFAULT_BOHRIUM_REGISTRY;
  const registryExists = fs.existsSync(registryPath);

  (registryExists ? it : it.skip)(
    'keeps OpenQC registry ids aligned with Bohrium backend ids',
    () => {
      const bohriumBackends = loadBohriumBackends(registryPath);
      const openqcIds = openqcServers.map(server => server.id).sort();
      const bohriumIds = bohriumBackends.map(backend => backend.id).sort();

      expect(openqcIds).toEqual(bohriumIds);
    }
  );

  (registryExists ? it : it.skip)(
    'keeps agent CLI names aligned between OpenQC readiness and Bohrium registry',
    () => {
      const bohriumById = new Map(
        loadBohriumBackends(registryPath).map(backend => [backend.id, backend])
      );

      for (const server of openqcServers) {
        const readiness = getLspDiagnosticReadiness(server.id);
        const backend = bohriumById.get(server.id);
        expect(backend).toBeDefined();
        expect(readiness?.agentCli).toBe(backend?.agent_cli);
      }
    }
  );

  it('exposes a local alignment script with JSON output', () => {
    const tempDir = fs.mkdtempSync(path.join(REPO_ROOT, '.tmp-bohrium-registry-'));
    const tempRegistry = path.join(tempDir, 'lsp_backends.yaml');
    const readinessById = new Map(
      openqcServers.map(server => [server.id, getLspDiagnosticReadiness(server.id)])
    );
    fs.writeFileSync(
      tempRegistry,
      JSON.stringify({
        backends: openqcServers.map(server => ({
          id: server.id,
          software: server.languageId,
          agent_cli: readinessById.get(server.id)?.agentCli,
        })),
      }),
      'utf8'
    );

    const output = execFileSync(
      process.execPath,
      [ALIGNMENT_SCRIPT, '--json', '--bohrium-registry', tempRegistry],
      {
        cwd: REPO_ROOT,
        encoding: 'utf8',
        stdio: ['ignore', 'pipe', 'pipe'],
      }
    );

    const report = JSON.parse(output) as {
      ok: boolean;
      surfaces: { openqc: { backendCount: number }; bohrium: { backendCount: number } };
      summary: { missingInBohrium: string[]; excessInBohrium: string[] };
    };

    expect(report.surfaces.openqc.backendCount).toBe(17);
    expect(report.surfaces.bohrium.backendCount).toBe(17);
    expect(report.summary.missingInBohrium).toEqual([]);
    expect(report.summary.excessInBohrium).toEqual([]);
    expect(report.ok).toBe(true);

    fs.rmSync(tempDir, { recursive: true, force: true });
  });
});
