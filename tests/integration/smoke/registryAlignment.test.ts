/**
 * Registry/Package/Docs Alignment Integration Tests
 *
 * Integration tests that verify alignment between the LSP registry,
 * package.json contributions, language configurations, grammars, and
 * documentation.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/164
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/165
 */

import * as fs from 'fs';
import * as path from 'path';
import { listBundledLspServers, getLspDiagnosticReadiness } from '../../../src/lsp/registry';
import { generateCompatibilityReport } from '../../../src/smoke/compatibilityReport';

const REPO_ROOT = path.resolve(__dirname, '../../../');

describe('Registry/Package/Docs Alignment Integration', () => {
  const servers = listBundledLspServers();
  const packageJson = JSON.parse(fs.readFileSync(path.join(REPO_ROOT, 'package.json'), 'utf8'));
  const languages = packageJson.contributes?.languages ?? [];
  const grammarDir = path.join(REPO_ROOT, 'syntaxes');
  const configDir = path.join(REPO_ROOT, 'language-configurations');
  const compatDoc = fs.readFileSync(path.join(REPO_ROOT, 'docs', 'LSP_COMPATIBILITY.md'), 'utf8');

  // -------------------------------------------------------------------------
  // Registry completeness
  // -------------------------------------------------------------------------

  test('every registry entry has a matching package.json language contribution', () => {
    for (const server of servers) {
      const lang = languages.find((l: Record<string, unknown>) => l.id === server.languageId);
      expect(lang).toBeDefined();

      // Verify extensions match
      const pkgExtensions = (lang?.extensions ?? []) as string[];
      for (const ext of server.fileExtensions) {
        expect(pkgExtensions).toContain(`.${ext}`);
      }

      // Verify filenames match
      const pkgFilenames = (lang?.filenames ?? []) as string[];
      for (const fn of server.fileNames) {
        expect(pkgFilenames).toContain(fn);
      }
    }
  });

  test('every package.json language has a registry entry', () => {
    for (const lang of languages) {
      const server = servers.find(s => s.languageId === lang.id);
      expect(server).toBeDefined();
    }
  });

  // -------------------------------------------------------------------------
  // Grammar files
  // -------------------------------------------------------------------------

  test('every registry entry has a grammar file', () => {
    const grammarFiles = fs.readdirSync(grammarDir);

    for (const server of servers) {
      const hasJsonGrammar = grammarFiles.includes(`${server.languageId}.tmLanguage.json`);
      const hasYamlGrammar = grammarFiles.includes(`${server.languageId}.tmLanguage.yaml`);
      const hasPlistGrammar = grammarFiles.includes(`${server.languageId}.tmLanguage`);

      expect(hasJsonGrammar || hasYamlGrammar || hasPlistGrammar).toBe(true);
    }
  });

  // -------------------------------------------------------------------------
  // Language configurations
  // -------------------------------------------------------------------------

  test('every registry entry has a language configuration file', () => {
    const configFiles = fs.readdirSync(configDir);

    for (const server of servers) {
      const hasConfig = configFiles.includes(`${server.languageId}.json`);
      expect(hasConfig).toBe(true);
    }
  });

  // -------------------------------------------------------------------------
  // Documentation alignment
  // -------------------------------------------------------------------------

  test('every registry entry is documented in LSP_COMPATIBILITY.md', () => {
    for (const server of servers) {
      expect(compatDoc).toContain(server.repository);
      expect(compatDoc).toContain(`\`${server.languageId}\``);
      expect(compatDoc).toContain(`\`${server.defaultBranch}\``);
    }
  });

  // -------------------------------------------------------------------------
  // Diagnostic readiness
  // -------------------------------------------------------------------------

  test('every registry entry has diagnostic readiness metadata', () => {
    for (const server of servers) {
      const readiness = getLspDiagnosticReadiness(server.id);
      expect(readiness).toBeDefined();
      expect(readiness!.diagnosticEngine).toBe('v1');
      expect(readiness!.richDiagnostics).toBe(true);
    }
  });

  // -------------------------------------------------------------------------
  // Configuration defaults alignment
  // -------------------------------------------------------------------------

  test('configuration defaults in package.json match registry', () => {
    const properties = packageJson.contributes?.configuration?.properties ?? {};

    for (const server of servers) {
      const enabledKey = `openqc.lsp.${server.languageId}.enabled`;
      const pathKey = `openqc.lsp.${server.languageId}.path`;
      const commandKey = `openqc.lsp.${server.languageId}.command`;

      expect(properties[enabledKey]).toBeDefined();
      expect(properties[enabledKey].default).toBe(server.enabled);

      expect(properties[pathKey]).toBeDefined();
      expect(properties[pathKey].default).toBe(server.executable);

      expect(properties[commandKey]).toBeDefined();
      expect(properties[commandKey].default).toBe(server.executable);
    }
  });

  // -------------------------------------------------------------------------
  // Full compatibility report
  // -------------------------------------------------------------------------

  test('full compatibility report passes for all servers', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    expect(report.totalServers).toBe(17);

    // All should pass core checks
    for (const entry of report.entries) {
      const registryCheck = entry.checks.find(c => c.name === 'registry-entry-exists');
      expect(registryCheck!.status).toBe('pass');

      const langCheck = entry.checks.find(c => c.name === 'package-json-language');
      expect(langCheck!.status).toBe('pass');

      const docsCheck = entry.checks.find(c => c.name === 'docs-alignment');
      expect(docsCheck!.status).toBe('pass');

      const diagCheck = entry.checks.find(c => c.name === 'diagnostic-readiness');
      expect(diagCheck!.status).toBe('pass');

      const launchCheck = entry.checks.find(c => c.name === 'local-launch-config');
      expect(launchCheck!.status).toBe('pass');
    }
  });
});
