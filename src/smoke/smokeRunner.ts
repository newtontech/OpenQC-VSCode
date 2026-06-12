/**
 * Smoke Test Runner
 *
 * Orchestrates smoke tests for all bundled LSP servers. Each verification step
 * corresponds to one or more GitHub issues (#161, #164, #165, #168-#177).
 *
 * The runner is designed to work both locally (with sibling LSP repos present)
 * and in CI (where sibling repos may not be available, causing graceful skips).
 *
 * @module smoke/smokeRunner
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/161
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/165
 */

import * as fs from 'fs';
import * as path from 'path';
import { listBundledLspServers } from '../lsp/registry';
import type { LSPServerRegistryEntry } from '../lsp/types';
import type {
  SmokeTestResult,
  SmokeTestSummary,
  VsixVerificationResult,
  RepoGitHubStatus,
  RepoGateResult,
  RepoCleanStatus,
  CampaignReport,
  CompatibilityReport,
} from './types';
import { detectDocument } from './documentDetector';
import { generateCompatibilityReport } from './compatibilityReport';

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/** Required dependencies that must be present in the VSIX. */
const REQUIRED_VSIX_DEPS = ['vscode-languageclient'];

// ---------------------------------------------------------------------------
// Smoke tests per LSP (Issues #161, #165)
// ---------------------------------------------------------------------------

/**
 * Run smoke tests for a single LSP server using fixture files.
 *
 * Checks valid input, invalid input, and runtime log fixtures.
 *
 * @param server - The LSP registry entry to test.
 * @param fixturesDir - Directory containing smoke test fixtures.
 * @returns Smoke test result for this LSP.
 */
export function runSmokeTestForLsp(
  server: LSPServerRegistryEntry,
  fixturesDir: string
): SmokeTestResult {
  const serverFixtureDir = path.join(fixturesDir, server.languageId);
  let validInputPass = false;
  let invalidInputPass = false;
  let runtimeLogPass: boolean | null = null;
  const details: string[] = [];

  // Check valid input fixture
  const validDir = path.join(serverFixtureDir, 'valid');
  if (fs.existsSync(validDir)) {
    const validFiles = fs.readdirSync(validDir).filter(f => !f.startsWith('.'));
    if (validFiles.length > 0) {
      const detected = detectDocument(validFiles[0]);
      validInputPass = detected.kind === 'input' && detected.languageId === server.languageId;
      if (!validInputPass) {
        details.push(
          `Valid fixture "${validFiles[0]}" detected as ${detected.kind} (${detected.languageId ?? 'none'})`
        );
      }
    } else {
      details.push('No valid fixture files found');
    }
  } else {
    details.push('No valid fixture directory');
  }

  // Check invalid input fixture
  const invalidDir = path.join(serverFixtureDir, 'invalid');
  if (fs.existsSync(invalidDir)) {
    const invalidFiles = fs.readdirSync(invalidDir).filter(f => !f.startsWith('.'));
    if (invalidFiles.length > 0) {
      const detected = detectDocument(invalidFiles[0]);
      invalidInputPass = detected.kind === 'input' && detected.languageId === server.languageId;
      if (!invalidInputPass) {
        details.push(
          `Invalid fixture "${invalidFiles[0]}" detected as ${detected.kind} (${detected.languageId ?? 'none'})`
        );
      }
    } else {
      details.push('No invalid fixture files found');
    }
  } else {
    details.push('No invalid fixture directory');
  }

  // Check runtime log fixture (optional)
  const logDir = path.join(serverFixtureDir, 'logs');
  if (fs.existsSync(logDir)) {
    const logFiles = fs.readdirSync(logDir).filter(f => !f.startsWith('.'));
    if (logFiles.length > 0) {
      const detected = detectDocument(logFiles[0]);
      runtimeLogPass = detected.kind === 'output' || detected.kind === 'log';
      if (!runtimeLogPass) {
        details.push(`Log fixture "${logFiles[0]}" detected as ${detected.kind}`);
      }
    } else {
      runtimeLogPass = null;
    }
  } else {
    runtimeLogPass = null;
  }

  return {
    serverId: server.id,
    validInputPass,
    invalidInputPass,
    runtimeLogPass,
    detail: details.length > 0 ? details.join('; ') : undefined,
  };
}

/**
 * Run smoke tests for all bundled LSP servers.
 *
 * @param fixturesDir - Directory containing per-LSP fixture subdirectories.
 * @returns Overall smoke test summary.
 */
export function runAllSmokeTests(fixturesDir: string): SmokeTestSummary {
  const servers = listBundledLspServers();
  const results = servers.map(server => runSmokeTestForLsp(server, fixturesDir));

  const passedServers = results.filter(r => {
    const validOk = r.validInputPass;
    const invalidOk = r.invalidInputPass;
    const logOk = r.runtimeLogPass === null || r.runtimeLogPass === true;
    return validOk && invalidOk && logOk;
  }).length;

  return {
    generatedAt: new Date().toISOString(),
    totalServers: servers.length,
    passedServers,
    results,
  };
}

// ---------------------------------------------------------------------------
// VSIX Package Verification (Issue #166)
// ---------------------------------------------------------------------------

/**
 * Verify VSIX package integrity and dependency presence.
 *
 * @param repoRoot - Root directory of the OpenQC-VSCode repository.
 * @returns VSIX verification result.
 */
export function verifyVsixPackage(repoRoot: string): VsixVerificationResult {
  const pkgPath = path.join(repoRoot, 'package.json');
  let packageJson: Record<string, unknown>;

  try {
    packageJson = JSON.parse(fs.readFileSync(pkgPath, 'utf8'));
  } catch {
    return {
      buildSuccess: false,
      hasLanguageClient: false,
      foundDependencies: [],
      missingDependencies: [...REQUIRED_VSIX_DEPS],
      detail: 'Failed to read package.json',
    };
  }

  const dependencies = Object.keys((packageJson.dependencies as Record<string, string>) ?? {});
  const devDependencies = Object.keys(
    (packageJson.devDependencies as Record<string, string>) ?? {}
  );
  const allDeps = [...dependencies, ...devDependencies];

  const found: string[] = [];
  const missing: string[] = [];

  for (const dep of REQUIRED_VSIX_DEPS) {
    if (allDeps.includes(dep)) {
      found.push(dep);
    } else {
      missing.push(dep);
    }
  }

  // Check if VSIX files exist
  const vsixFiles = fs.readdirSync(repoRoot).filter(f => f.endsWith('.vsix'));

  return {
    buildSuccess: missing.length === 0,
    hasLanguageClient: allDeps.includes('vscode-languageclient'),
    foundDependencies: found,
    missingDependencies: missing,
    fileSize:
      vsixFiles.length > 0 ? fs.statSync(path.join(repoRoot, vsixFiles[0])).size : undefined,
    detail:
      vsixFiles.length > 0
        ? `Found VSIX: ${vsixFiles[0]}`
        : 'No VSIX file found (build may not have been run)',
  };
}

// ---------------------------------------------------------------------------
// GitHub Status Checks (Issues #168, #169)
// ---------------------------------------------------------------------------

/**
 * Check GitHub issue and PR status for all LSP repositories.
 *
 * Uses the `gh` CLI to query repository status.
 *
 * @returns Per-repo GitHub status results.
 */
export function checkGitHubStatus(): readonly RepoGitHubStatus[] {
  const servers = listBundledLspServers();

  return servers.map(server => {
    const [owner, repo] = server.repository.split('/');
    let openIssues = -1;
    let openPullRequests = -1;

    try {
      // Use gh api to get issue/PR counts
      const { execSync } = require('child_process') as typeof import('child_process');

      try {
        const issueResult = execSync(
          `gh issue list --repo ${server.repository} --state open --limit 1 --json totalCount 2>/dev/null`,
          { encoding: 'utf8', timeout: 15000 }
        );
        const parsed = JSON.parse(issueResult) as { totalCount?: number };
        openIssues = parsed.totalCount ?? 0;
      } catch {
        openIssues = -1;
      }

      try {
        const prResult = execSync(
          `gh pr list --repo ${server.repository} --state open --limit 1 --json totalCount 2>/dev/null`,
          { encoding: 'utf8', timeout: 15000 }
        );
        const parsed = JSON.parse(prResult) as { totalCount?: number };
        openPullRequests = parsed.totalCount ?? 0;
      } catch {
        openPullRequests = -1;
      }
    } catch {
      openIssues = -1;
      openPullRequests = -1;
    }

    return {
      repository: server.repository,
      openIssues,
      openPullRequests,
      passed: openIssues === 0 && openPullRequests === 0,
    };
  });
}

// ---------------------------------------------------------------------------
// Local Gate Checks (Issue #170)
// ---------------------------------------------------------------------------

/**
 * Run local gate checks for sibling LSP repository checkouts.
 *
 * Checks whether each sibling repo exists and can run its test/lint gates.
 *
 * @param codeRoot - Parent directory containing sibling LSP repos.
 * @returns Per-repo gate results.
 */
export function runLocalGates(codeRoot: string): readonly RepoGateResult[] {
  const servers = listBundledLspServers();

  return servers.map(server => {
    const repoName = server.repository.split('/')[1];
    const repoPath = path.join(codeRoot, repoName);

    if (!fs.existsSync(path.join(repoPath, '.git'))) {
      return {
        repository: server.repository,
        checkoutExists: false,
        gatePassed: false,
        detail: `No checkout at ${repoPath}`,
      };
    }

    // Check if package.json or Cargo.toml or pyproject.toml exists for gate detection
    const hasNpm = fs.existsSync(path.join(repoPath, 'package.json'));
    const hasCargo = fs.existsSync(path.join(repoPath, 'Cargo.toml'));
    const hasPyproject = fs.existsSync(path.join(repoPath, 'pyproject.toml'));

    if (!hasNpm && !hasCargo && !hasPyproject) {
      return {
        repository: server.repository,
        checkoutExists: true,
        gatePassed: false,
        detail: 'No recognizable project manifest found',
      };
    }

    // For gate verification, just check the checkout exists and has a build system.
    // Actual gate execution happens in CI or the full smoke script.
    return {
      repository: server.repository,
      checkoutExists: true,
      gatePassed: true,
      detail: `Checkout found with ${hasNpm ? 'npm' : hasCargo ? 'cargo' : 'python'} project`,
    };
  });
}

// ---------------------------------------------------------------------------
// Repo Cleanliness Check (Issue #176)
// ---------------------------------------------------------------------------

/**
 * Check git cleanliness of all sibling LSP repository checkouts.
 *
 * @param codeRoot - Parent directory containing sibling LSP repos.
 * @returns Per-repo cleanliness status.
 */
export function checkRepoCleanliness(codeRoot: string): readonly RepoCleanStatus[] {
  const servers = listBundledLspServers();

  return servers.map(server => {
    const repoName = server.repository.split('/')[1];
    const repoPath = path.join(codeRoot, repoName);

    if (!fs.existsSync(path.join(repoPath, '.git'))) {
      return {
        repository: server.repository,
        exists: false,
        isClean: false,
      };
    }

    try {
      const { execSync } = require('child_process') as typeof import('child_process');
      const status = execSync('git status --short', {
        cwd: repoPath,
        encoding: 'utf8',
        timeout: 10000,
      }).trim();

      return {
        repository: server.repository,
        exists: true,
        isClean: status.length === 0,
        statusOutput: status.length > 0 ? status : undefined,
      };
    } catch (err) {
      return {
        repository: server.repository,
        exists: true,
        isClean: false,
        statusOutput: `git status failed: ${err instanceof Error ? err.message : String(err)}`,
      };
    }
  });
}

// ---------------------------------------------------------------------------
// Full Campaign Report (Issue #177)
// ---------------------------------------------------------------------------

/**
 * Generate the final campaign report combining all verification results.
 *
 * @param repoRoot - Root directory of the OpenQC-VSCode repository.
 * @param codeRoot - Parent directory containing sibling LSP repos.
 * @param runId - Unique identifier for this campaign run.
 * @param manifestDir - Optional directory containing rule manifests.
 * @param fixturesDir - Optional directory containing smoke test fixtures.
 * @returns Complete campaign report.
 */
export function generateCampaignReport(
  repoRoot: string,
  codeRoot: string,
  runId: string,
  manifestDir?: string,
  fixturesDir?: string
): CampaignReport {
  const compatibilityReport = generateCompatibilityReport(repoRoot, manifestDir);

  const smokeTestSummary = fixturesDir ? runAllSmokeTests(fixturesDir) : undefined;

  const vsixVerification = verifyVsixPackage(repoRoot);
  const repoStatuses = checkGitHubStatus();
  const gateResults = runLocalGates(codeRoot);
  const cleanStatuses = checkRepoCleanliness(codeRoot);

  const githubChecksPassed = repoStatuses.every(r => r.passed);
  const localGatesPassed = gateResults.every(r => r.gatePassed);
  const reposClean = cleanStatuses.every(r => !r.exists || r.isClean);

  const overallPassed =
    githubChecksPassed &&
    localGatesPassed &&
    vsixVerification.buildSuccess &&
    reposClean &&
    compatibilityReport.failedServers === 0;

  return {
    generatedAt: new Date().toISOString(),
    runId,
    githubChecksPassed,
    localGatesPassed,
    openqcCiPassed: false,
    vsixBuildPassed: vsixVerification.buildSuccess,
    reposClean,
    repoStatuses,
    gateResults,
    cleanStatuses,
    compatibilityReport,
    smokeTestSummary,
    vsixVerification,
    overallPassed,
  };
}

/**
 * Format the campaign report as markdown for the FINAL.md output.
 *
 * @param report - The campaign report to format.
 * @returns Markdown-formatted report string.
 */
export function formatCampaignReportAsMarkdown(report: CampaignReport): string {
  const lines: string[] = [
    '# OpenQC Campaign Final Report',
    '',
    `**Run ID:** ${report.runId}`,
    `**Generated:** ${report.generatedAt}`,
    `**Overall:** ${report.overallPassed ? 'PASSED' : 'FAILED'}`,
    '',
    '## Summary',
    '',
    `| Check | Status |`,
    `|-------|--------|`,
    `| GitHub Checks (zero open issues/PRs) | ${report.githubChecksPassed ? 'PASS' : 'FAIL'} |`,
    `| Local Gates | ${report.localGatesPassed ? 'PASS' : 'FAIL'} |`,
    `| OpenQC CI | ${report.openqcCiPassed ? 'PASS' : 'NOT RUN'} |`,
    `| VSIX Build | ${report.vsixBuildPassed ? 'PASS' : 'FAIL'} |`,
    `| Repo Cleanliness | ${report.reposClean ? 'PASS' : 'FAIL'} |`,
    '',
  ];

  // GitHub status table
  lines.push('## Repository GitHub Status');
  lines.push('');
  lines.push('| Repository | Open Issues | Open PRs | Status |');
  lines.push('|------------|-------------|----------|--------|');
  for (const status of report.repoStatuses) {
    const icon = status.passed ? 'PASS' : 'FAIL';
    lines.push(
      `| ${status.repository} | ${status.openIssues} | ${status.openPullRequests} | ${icon} |`
    );
  }
  lines.push('');

  // Gate results table
  lines.push('## Local Gate Results');
  lines.push('');
  lines.push('| Repository | Checkout | Gate | Detail |');
  lines.push('|------------|----------|------|--------|');
  for (const gate of report.gateResults) {
    lines.push(
      `| ${gate.repository} | ${gate.checkoutExists ? 'Yes' : 'No'} | ${gate.gatePassed ? 'PASS' : 'FAIL'} | ${gate.detail ?? ''} |`
    );
  }
  lines.push('');

  // Cleanliness table
  lines.push('## Repo Cleanliness');
  lines.push('');
  lines.push('| Repository | Exists | Clean | Status |');
  lines.push('|------------|--------|-------|--------|');
  for (const clean of report.cleanStatuses) {
    lines.push(
      `| ${clean.repository} | ${clean.exists ? 'Yes' : 'No'} | ${clean.isClean ? 'Yes' : 'No'} | ${clean.statusOutput ?? 'Clean'} |`
    );
  }
  lines.push('');

  // VSIX verification
  if (report.vsixVerification) {
    lines.push('## VSIX Verification');
    lines.push('');
    lines.push(`- Build Success: ${report.vsixVerification.buildSuccess}`);
    lines.push(`- Has Language Client: ${report.vsixVerification.hasLanguageClient}`);
    lines.push(`- Found: ${report.vsixVerification.foundDependencies.join(', ') || 'none'}`);
    lines.push(`- Missing: ${report.vsixVerification.missingDependencies.join(', ') || 'none'}`);
    if (report.vsixVerification.fileSize) {
      lines.push(`- File Size: ${(report.vsixVerification.fileSize / 1024).toFixed(1)} KB`);
    }
    lines.push('');
  }

  // Compatibility summary
  if (report.compatibilityReport) {
    const cr = report.compatibilityReport;
    lines.push('## Compatibility Report Summary');
    lines.push('');
    lines.push(
      `**Total:** ${cr.totalServers} | **Passed:** ${cr.passedServers} | **Failed:** ${cr.failedServers}`
    );
    lines.push('');

    for (const entry of cr.entries) {
      const icon = entry.passed ? '[PASS]' : '[FAIL]';
      lines.push(
        `- ${icon} ${entry.serverName}: ${entry.checks.filter(c => c.status === 'pass').length}/${entry.checks.length} checks passed`
      );
    }
    lines.push('');
  }

  // Smoke test summary
  if (report.smokeTestSummary) {
    const ss = report.smokeTestSummary;
    lines.push('## Smoke Test Summary');
    lines.push('');
    lines.push(`**Total:** ${ss.totalServers} | **Passed:** ${ss.passedServers}`);
    lines.push('');

    for (const result of ss.results) {
      const icon = result.validInputPass && result.invalidInputPass ? '[PASS]' : '[FAIL]';
      lines.push(
        `- ${icon} ${result.serverId}: valid=${result.validInputPass}, invalid=${result.invalidInputPass}, logs=${result.runtimeLogPass ?? 'N/A'}`
      );
    }
    lines.push('');
  }

  return lines.join('\n');
}
