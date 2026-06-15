/**
 * Compatibility Report Tests
 *
 * Tests for the compatibility report generator.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/163
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/167
 */

import * as path from 'path';
import {
  generateCompatibilityReport,
  formatReportAsMarkdown,
} from '../../../src/smoke/compatibilityReport';
import { getBundledLspServerCount } from '../../../src/lsp/registry';

const REPO_ROOT = path.resolve(__dirname, '../../../');
const MANIFESTS_DIR = path.join(REPO_ROOT, 'tests/fixtures/smoke/manifests');

// ---------------------------------------------------------------------------
// generateCompatibilityReport
// ---------------------------------------------------------------------------

describe('generateCompatibilityReport', () => {
  it('generates a report covering every bundled LSP server', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    expect(report.totalServers).toBe(getBundledLspServerCount());
    expect(report.entries).toHaveLength(getBundledLspServerCount());
  });

  it('includes timestamps', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    expect(report.generatedAt).toBeTruthy();
    expect(new Date(report.generatedAt).getTime()).not.toBeNaN();
  });

  it('counts passed and failed servers', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    expect(report.passedServers + report.failedServers).toBe(report.totalServers);
  });

  it('checks registry entry for each server', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    for (const entry of report.entries) {
      const registryCheck = entry.checks.find(c => c.name === 'registry-entry-exists');
      expect(registryCheck).toBeDefined();
      expect(registryCheck!.status).toBe('pass');
    }
  });

  it('checks package.json language contribution for each server', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    for (const entry of report.entries) {
      const langCheck = entry.checks.find(c => c.name === 'package-json-language');
      expect(langCheck).toBeDefined();
      expect(langCheck!.status).toBe('pass');
    }
  });

  it('checks package.json configuration for each server', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    for (const entry of report.entries) {
      const configCheck = entry.checks.find(c => c.name === 'package-json-configuration');
      expect(configCheck).toBeDefined();
      // Should pass or warn (some may be missing optional settings)
      expect(['pass', 'warn']).toContain(configCheck!.status);
    }
  });

  it('checks docs alignment for each server', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    for (const entry of report.entries) {
      const docsCheck = entry.checks.find(c => c.name === 'docs-alignment');
      expect(docsCheck).toBeDefined();
      expect(docsCheck!.status).toBe('pass');
    }
  });

  it('checks diagnostic readiness for each server', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    for (const entry of report.entries) {
      const diagCheck = entry.checks.find(c => c.name === 'diagnostic-readiness');
      expect(diagCheck).toBeDefined();
      expect(diagCheck!.status).toBe('pass');
    }
  });

  it('checks local launch config for each server', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    for (const entry of report.entries) {
      const launchCheck = entry.checks.find(c => c.name === 'local-launch-config');
      expect(launchCheck).toBeDefined();
      expect(launchCheck!.status).toBe('pass');
    }
  });

  it('checks capability manifests for sibling LSP repositories that expose them', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    for (const entry of report.entries) {
      const manifestCheck = entry.checks.find(c => c.name === 'lsp-capability-manifest');
      expect(manifestCheck).toBeDefined();
      if (manifestCheck!.status === 'pass') {
        expect(entry.manifest).toBeDefined();
        expect(entry.manifest!.agentCli).toBeTruthy();
        expect(entry.manifest!.capabilities).toContain('diagnostics');
      }
      expect(entry.runtimeCommand).toContain(entry.languageId === 'cp2k' ? 'cp2k' : '');
    }
  });

  it('includes runtime and manifest fields in each report entry', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    const dpgen = report.entries.find(e => e.serverId === 'dpgen-lsp');
    expect(dpgen).toBeDefined();
    expect(dpgen!.runtimeCommand).toBe('dpgen-lsp --stdio');
    expect(dpgen!.repoPath).toContain('dpgen-lsp');
    expect(dpgen!.manifest?.blockingPolicy).toBe('blocking');
    expect(dpgen!.manifest?.missingCapabilities).toEqual([]);
  });

  it('checks rule manifests when manifestDir is provided', () => {
    const report = generateCompatibilityReport(REPO_ROOT, MANIFESTS_DIR);
    const gaussian = report.entries.find(e => e.serverId === 'gaussian-lsp');
    expect(gaussian).toBeDefined();
    const manifestCheck = gaussian!.checks.find(c => c.name === 'rule-manifest');
    expect(manifestCheck).toBeDefined();
    expect(manifestCheck!.status).toBe('pass');
    expect(manifestCheck!.detail).toContain('4 rules');
  });

  it('skips rule manifest check when no manifestDir provided', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    for (const entry of report.entries) {
      const manifestCheck = entry.checks.find(c => c.name === 'rule-manifest');
      expect(manifestCheck).toBeDefined();
      expect(manifestCheck!.status).toBe('skip');
    }
  });

  it('handles invalid repoRoot gracefully', () => {
    const report = generateCompatibilityReport('/nonexistent/path');
    expect(report.totalServers).toBe(getBundledLspServerCount());
    // Should still generate the report but with failures
    for (const entry of report.entries) {
      const docsCheck = entry.checks.find(c => c.name === 'docs-alignment');
      expect(docsCheck!.status).toBe('fail');
    }
  });
});

// ---------------------------------------------------------------------------
// formatReportAsMarkdown
// ---------------------------------------------------------------------------

describe('formatReportAsMarkdown', () => {
  it('produces markdown output', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    const md = formatReportAsMarkdown(report);
    expect(md).toContain('# OpenQC LSP Compatibility Report');
    expect(md).toContain('## Per-Server Results');
  });

  it('includes server names in the markdown', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    const md = formatReportAsMarkdown(report);
    expect(md).toContain('Gaussian');
    expect(md).toContain('VASP');
  });

  it('includes check details', () => {
    const report = generateCompatibilityReport(REPO_ROOT);
    const md = formatReportAsMarkdown(report);
    expect(md).toContain('registry-entry-exists');
    expect(md).toContain('package-json-language');
    expect(md).toContain('lsp-capability-manifest');
    expect(md).toContain('Runtime command');
  });
});
