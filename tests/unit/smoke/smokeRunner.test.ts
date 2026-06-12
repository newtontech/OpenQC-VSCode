/**
 * Smoke Runner Tests
 *
 * Tests for the smoke test runner, VSIX verification, and campaign report.
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/161
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/164
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/165
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/166
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/168
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/176
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/177
 */

import * as path from 'path';
import {
  runSmokeTestForLsp,
  runAllSmokeTests,
  verifyVsixPackage,
  runLocalGates,
  checkRepoCleanliness,
  generateCampaignReport,
  formatCampaignReportAsMarkdown,
} from '../../../src/smoke/smokeRunner';
import { getLspServerByLanguageId } from '../../../src/lsp/registry';

const REPO_ROOT = path.resolve(__dirname, '../../../');
const CODE_ROOT = path.resolve(REPO_ROOT, '..');

// ---------------------------------------------------------------------------
// runSmokeTestForLsp
// ---------------------------------------------------------------------------

describe('runSmokeTestForLsp', () => {
  it('returns a result for a valid server with no fixtures', () => {
    const server = getLspServerByLanguageId('gaussian')!;
    const result = runSmokeTestForLsp(server, '/nonexistent/fixtures');
    expect(result.serverId).toBe('gaussian-lsp');
    expect(result.validInputPass).toBe(false);
    expect(result.invalidInputPass).toBe(false);
    expect(result.runtimeLogPass).toBe(null);
    expect(result.detail).toContain('No valid fixture');
  });
});

// ---------------------------------------------------------------------------
// runAllSmokeTests
// ---------------------------------------------------------------------------

describe('runAllSmokeTests', () => {
  it('runs smoke tests for all 16 servers', () => {
    const summary = runAllSmokeTests('/nonexistent/fixtures');
    expect(summary.totalServers).toBe(16);
    expect(summary.results).toHaveLength(16);
    expect(summary.generatedAt).toBeTruthy();
  });
});

// ---------------------------------------------------------------------------
// verifyVsixPackage
// ---------------------------------------------------------------------------

describe('verifyVsixPackage', () => {
  it('checks for required dependencies in package.json', () => {
    const result = verifyVsixPackage(REPO_ROOT);
    expect(result.hasLanguageClient).toBe(true);
    expect(result.foundDependencies).toContain('vscode-languageclient');
    expect(result.missingDependencies).toEqual([]);
    expect(result.buildSuccess).toBe(true);
  });

  it('reports failure for invalid repo root', () => {
    const result = verifyVsixPackage('/nonexistent/path');
    expect(result.buildSuccess).toBe(false);
    expect(result.hasLanguageClient).toBe(false);
  });
});

// ---------------------------------------------------------------------------
// runLocalGates
// ---------------------------------------------------------------------------

describe('runLocalGates', () => {
  it('checks all 16 repos', () => {
    const results = runLocalGates(CODE_ROOT);
    expect(results).toHaveLength(16);
  });

  it('reports missing checkouts gracefully', () => {
    const results = runLocalGates('/nonexistent/code');
    for (const result of results) {
      expect(result.checkoutExists).toBe(false);
      expect(result.gatePassed).toBe(false);
    }
  });
});

// ---------------------------------------------------------------------------
// checkRepoCleanliness
// ---------------------------------------------------------------------------

describe('checkRepoCleanliness', () => {
  it('checks all 16 repos', () => {
    const statuses = checkRepoCleanliness(CODE_ROOT);
    expect(statuses).toHaveLength(16);
  });

  it('reports non-existent repos', () => {
    const statuses = checkRepoCleanliness('/nonexistent/code');
    for (const status of statuses) {
      expect(status.exists).toBe(false);
    }
  });
});

// ---------------------------------------------------------------------------
// generateCampaignReport
// ---------------------------------------------------------------------------

describe('generateCampaignReport', () => {
  it('generates a full campaign report', () => {
    const report = generateCampaignReport(
      REPO_ROOT,
      CODE_ROOT,
      'test-run-001',
      path.join(REPO_ROOT, 'tests/fixtures/smoke/manifests')
    );

    expect(report.runId).toBe('test-run-001');
    expect(report.generatedAt).toBeTruthy();
    expect(report.repoStatuses).toHaveLength(16);
    expect(report.gateResults).toHaveLength(16);
    expect(report.cleanStatuses).toHaveLength(16);
    expect(report.compatibilityReport).toBeDefined();
    if (report.compatibilityReport) {
      expect(report.compatibilityReport.totalServers).toBe(16);
    }
  });

  it('includes VSIX verification in report', () => {
    const report = generateCampaignReport(REPO_ROOT, CODE_ROOT, 'test-run-002');
    expect(report.vsixVerification).toBeDefined();
    expect(report.vsixVerification!.hasLanguageClient).toBe(true);
  });

  it('works without manifestDir or fixturesDir', () => {
    const report = generateCampaignReport(REPO_ROOT, CODE_ROOT, 'test-run-003');
    expect(report.compatibilityReport).toBeDefined();
    expect(report.smokeTestSummary).toBeUndefined();
  });
});

// ---------------------------------------------------------------------------
// formatCampaignReportAsMarkdown
// ---------------------------------------------------------------------------

describe('formatCampaignReportAsMarkdown', () => {
  it('produces markdown with all sections', () => {
    const report = generateCampaignReport(
      REPO_ROOT,
      CODE_ROOT,
      'test-md',
      path.join(REPO_ROOT, 'tests/fixtures/smoke/manifests'),
      '/nonexistent/fixtures'
    );
    const md = formatCampaignReportAsMarkdown(report);

    expect(md).toContain('# OpenQC Campaign Final Report');
    expect(md).toContain('**Run ID:** test-md');
    expect(md).toContain('## Repository GitHub Status');
    expect(md).toContain('## Local Gate Results');
    expect(md).toContain('## Repo Cleanliness');
    expect(md).toContain('## VSIX Verification');
    expect(md).toContain('## Compatibility Report Summary');
    expect(md).toContain('## Smoke Test Summary');
  });

  it('includes server names in the report', () => {
    const report = generateCampaignReport(REPO_ROOT, CODE_ROOT, 'test-names');
    const md = formatCampaignReportAsMarkdown(report);
    expect(md).toContain('Gaussian');
    expect(md).toContain('VASP');
  });
});
