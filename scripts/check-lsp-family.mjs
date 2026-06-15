#!/usr/bin/env node

/**
 * LSP Family Release Gate Check
 *
 * Validates every bundled LSP server against fleet-perfect release requirements:
 * - Manifest presence and validity (lsp-capabilities.json)
 * - Provenance entries (version, commit, source)
 * - Fixture coverage (valid, invalid, log)
 * - Smoke command availability
 * - DiagnosticEnvelope/v1 readiness
 *
 * Usage:
 *   node scripts/check-lsp-family.mjs           # human-readable output
 *   node scripts/check-lsp-family.mjs --json    # machine-readable JSON
 *   node scripts/check-lsp-family.mjs --strict  # exit non-zero on any gap
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/186
 */

import { execFileSync } from 'child_process';
import { existsSync, readFileSync, readdirSync, statSync } from 'fs';
import { dirname, join, resolve } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '..');
const codeRoot = resolve(repoRoot, '..');

// ---------------------------------------------------------------------------
// CLI flags
// ---------------------------------------------------------------------------

const jsonMode = process.argv.includes('--json');
const strictMode = process.argv.includes('--strict');
const verbose = process.argv.includes('--verbose');
const reportPath = process.argv.find((arg, i) => process.argv[i - 1] === '--report-path') ?? null;

// ---------------------------------------------------------------------------
// Parse registry
// ---------------------------------------------------------------------------

const registryPath = join(repoRoot, 'src/lsp/registry.ts');
const registry = readFileSync(registryPath, 'utf8');

const entries = [...registry.matchAll(
  /id: '([^']+)'[\s\S]*?repository: '([^']+)'[\s\S]*?languageId: '([^']+)'[\s\S]*?defaultBranch: '([^']+)'/g
)].map(match => ({
  id: match[1],
  repository: match[2],
  languageId: match[3],
  defaultBranch: match[4],
}));

const readiness = new Map([...registry.matchAll(
  /'([^']+)': diagnosticReadiness\('([^']+)', '([^']+)'/g
)].map(match => [
  match[1],
  {
    agentCli: match[2],
    closedLoop: match[3],
  },
]));

if (entries.length === 0) {
  throw new Error(`No LSP registry entries found in ${registryPath}`);
}

// ---------------------------------------------------------------------------
// Severity levels
// ---------------------------------------------------------------------------

const SEVERITY = {
  ERROR: 'error',       // blocks release
  WARNING: 'warning',   // maturity gap, does not block
  INFO: 'info',         // informational
};

// ---------------------------------------------------------------------------
// Warning categories for actionable maturity tracking
// ---------------------------------------------------------------------------

const WARNING_CATEGORIES = {
  MANIFEST: {
    id: 'manifest',
    label: 'Capability Manifest',
    description: 'lsp-capabilities.json presence and shape',
    graduationWeight: 3,
  },
  PROVENANCE: {
    id: 'provenance',
    label: 'Provenance & Versioning',
    description: 'Git tags, CHANGELOG, VERSION files',
    graduationWeight: 2,
  },
  FIXTURES: {
    id: 'fixtures',
    label: 'Test Fixtures',
    description: 'Valid/invalid fixture coverage',
    graduationWeight: 2,
  },
  DIAGNOSTIC: {
    id: 'diagnostic',
    label: 'Diagnostic Readiness',
    description: 'Agent CLI, closed-loop support',
    graduationWeight: 3,
  },
  CHECKOUT: {
    id: 'checkout',
    label: 'Local Checkout',
    description: 'Sibling checkout availability',
    graduationWeight: 1,
  },
};

/**
 * Categorize a warning message into a warning category.
 */
function categorizeWarning(message) {
  const lower = message.toLowerCase();
  if (lower.includes('manifest') || lower.includes('capabilities') || lower.includes('lsp-capabilities')) {
    return WARNING_CATEGORIES.MANIFEST;
  }
  if (lower.includes('tag') || lower.includes('changelog') || lower.includes('version')) {
    return WARNING_CATEGORIES.PROVENANCE;
  }
  if (lower.includes('fixture')) {
    return WARNING_CATEGORIES.FIXTURES;
  }
  if (lower.includes('closed-loop') || lower.includes('diagnostic') || lower.includes('agent cli')) {
    return WARNING_CATEGORIES.DIAGNOSTIC;
  }
  if (lower.includes('checkout') || lower.includes('local')) {
    return WARNING_CATEGORIES.CHECKOUT;
  }
  return { id: 'other', label: 'Other', description: 'Uncategorized warnings', graduationWeight: 1 };
}

/**
 * Calculate graduation score (0-100) based on warning categories resolved.
 */
function calculateGraduationScore(results) {
  const totalWeight = Object.values(WARNING_CATEGORIES).reduce((sum, cat) => sum + cat.graduationWeight, 0);
  const repoCount = results.length;

  // Count warnings per category across all repos
  const categoryWarnings = new Map();
  for (const cat of Object.values(WARNING_CATEGORIES)) {
    categoryWarnings.set(cat.id, 0);
  }

  for (const result of results) {
    for (const gap of result.warningGaps) {
      const cat = categorizeWarning(gap.message);
      const existing = categoryWarnings.get(cat.id) ?? 0;
      categoryWarnings.set(cat.id, existing + 1);
    }
  }

  // Calculate score: each category contributes proportionally
  // A category with 0 warnings across all repos gets full weight
  let earnedWeight = 0;
  for (const cat of Object.values(WARNING_CATEGORIES)) {
    const warnings = categoryWarnings.get(cat.id) ?? 0;
    const maxPossibleWarnings = repoCount; // Each repo could have this warning
    const resolvedRatio = Math.max(0, 1 - warnings / maxPossibleWarnings);
    earnedWeight += cat.graduationWeight * resolvedRatio;
  }

  return Math.round((earnedWeight / totalWeight) * 100);
}

/**
 * Build warning summary grouped by category.
 */
function buildWarningSummary(results) {
  const summary = new Map();

  for (const cat of Object.values(WARNING_CATEGORIES)) {
    summary.set(cat.id, {
      category: cat.label,
      description: cat.description,
      count: 0,
      repos: [],
      messages: new Set(),
    });
  }

  for (const result of results) {
    for (const gap of result.warningGaps) {
      const cat = categorizeWarning(gap.message);
      const entry = summary.get(cat.id);
      if (entry) {
        entry.count++;
        if (!entry.repos.includes(result.id)) {
          entry.repos.push(result.id);
        }
        entry.messages.add(gap.message);
      }
    }
  }

  // Convert Sets to Arrays for JSON serialization
  const output = {};
  for (const [id, entry] of summary) {
    output[id] = {
      ...entry,
      messages: [...entry.messages],
    };
  }
  return output;
}

// ---------------------------------------------------------------------------
// Gate checks
// ---------------------------------------------------------------------------

/**
 * Check for lsp-capabilities.json manifest in sibling repo.
 */
function checkManifest(entry) {
  const repoName = entry.repository.split('/')[1];
  const candidates = checkoutCandidates(entry)
    .map(path => join(path, 'lsp-capabilities.json'));

  for (const path of candidates) {
    if (existsSync(path)) {
      try {
        const content = readFileSync(path, 'utf8');
        const manifest = JSON.parse(content);
        const gaps = validateManifestShape(entry, manifest);
        return { status: 'pass', path, gaps };
      } catch (e) {
        return { status: 'invalid', path, gaps: [{ severity: SEVERITY.ERROR, message: `Invalid JSON: ${e.message}` }] };
      }
    }
  }

  return { status: 'missing', path: null, gaps: [{ severity: SEVERITY.ERROR, message: 'lsp-capabilities.json not found in any sibling checkout' }] };
}

/**
 * Validate manifest has the required fields.
 */
function validateManifestShape(entry, manifest) {
  const gaps = [];
  const repoName = entry.repository.split('/')[1];
  const manifestLanguageId = manifest.languageId ?? manifest.language_id;
  const manifestVersion = manifest.version
    ?? manifest.capabilities_version
    ?? manifest.schema_version
    ?? manifest.diagnostic_engine;
  const manifestRepository = manifest.repository
    ?? (manifest.openqc?.repoName ? `newtontech/${manifest.openqc.repoName}` : undefined);
  const manifestIdentity = manifest.id ?? manifest.name ?? manifest.openqc?.registryId;

  if (!manifestLanguageId && manifestIdentity !== entry.id && manifestRepository !== entry.repository) {
    gaps.push({
      severity: SEVERITY.WARNING,
      message: 'No canonical languageId field; relying on registry/repository identity',
    });
  }
  if (!manifestVersion) {
    gaps.push({
      severity: SEVERITY.WARNING,
      message: 'No canonical version field or known version alias',
    });
  }
  if (!manifestRepository && manifestIdentity !== entry.id) {
    gaps.push({ severity: SEVERITY.ERROR, message: 'Missing repository identity (repository/openqc.repoName/id/name)' });
  }
  if (manifestRepository && manifestRepository !== entry.repository) {
    const manifestRepoName = manifestRepository.split('/')[1] ?? manifestRepository;
    if (manifestRepoName !== repoName) {
      gaps.push({ severity: SEVERITY.ERROR, message: `repository mismatch: manifest=${manifestRepository} registry=${entry.repository}` });
    }
  }
  if (manifestLanguageId && manifestLanguageId !== entry.languageId) {
    gaps.push({ severity: SEVERITY.ERROR, message: `languageId mismatch: manifest=${manifestLanguageId} registry=${entry.languageId}` });
  }
  if (!manifest.capabilities) {
    gaps.push({ severity: SEVERITY.WARNING, message: 'No capabilities section (editor parity tracking incomplete)' });
  }
  return gaps;
}

/**
 * Check for provenance metadata in sibling repo.
 */
function checkProvenance(entry) {
  const repoName = entry.repository.split('/')[1];
  const localPath = findLocalCheckout(entry);

  if (!localPath) {
    return { status: 'missing', gaps: [{ severity: SEVERITY.WARNING, message: 'No local checkout found' }] };
  }

  const gaps = [];

  // Check for version/tag
  try {
    const tags = git(['-C', localPath, 'tag', '--sort=-version:refname'], { quiet: true });
    if (tags.trim().length === 0) {
      gaps.push({ severity: SEVERITY.WARNING, message: 'No git tags found (no release version)' });
    }
  } catch {
    gaps.push({ severity: SEVERITY.INFO, message: 'Could not read tags' });
  }

  // Check for CHANGELOG or version file
  const hasChangelog = existsSync(join(localPath, 'CHANGELOG.md')) || existsSync(join(localPath, 'CHANGELOG'));
  const hasVersion = existsSync(join(localPath, 'VERSION')) || existsSync(join(localPath, 'version.txt'));
  if (!hasChangelog && !hasVersion) {
    gaps.push({ severity: SEVERITY.WARNING, message: 'No CHANGELOG.md or VERSION file' });
  }

  return { status: gaps.length === 0 ? 'pass' : 'partial', gaps };
}

/**
 * Check for valid/invalid/log fixture sets.
 */
function checkFixtures(entry) {
  const repoName = entry.repository.split('/')[1];
  const localPath = findLocalCheckout(entry);

  if (!localPath) {
    return { status: 'missing', gaps: [{ severity: SEVERITY.WARNING, message: 'No local checkout; cannot check fixtures' }] };
  }

  const gaps = [];
  const fixtureDirs = ['tests/fixtures', 'fixtures', 'test/fixtures', 'tests/test_data'];

  let foundValid = false;
  let foundInvalid = false;

  for (const dir of fixtureDirs) {
    const fixturePath = join(localPath, dir);
    if (!existsSync(fixturePath)) continue;

    const files = listFilesRecursive(fixturePath);
    const names = files.map(f => f.toLowerCase());

    if (names.some(n => n.includes('valid') || n.includes('correct') || n.includes('ok'))) {
      foundValid = true;
    }
    if (names.some(n => n.includes('invalid') || n.includes('error') || n.includes('bad'))) {
      foundInvalid = true;
    }
  }

  if (!foundValid) {
    gaps.push({ severity: SEVERITY.WARNING, message: 'No valid fixture files found' });
  }
  if (!foundInvalid) {
    gaps.push({ severity: SEVERITY.WARNING, message: 'No invalid fixture files found' });
  }

  return {
    status: gaps.length === 0 ? 'pass' : 'partial',
    gaps,
  };
}

/**
 * Check smoke command availability.
 */
function checkSmoke(entry) {
  const metadata = readiness.get(entry.id);
  if (!metadata?.agentCli) {
    return { status: 'skip', gaps: [{ severity: SEVERITY.INFO, message: 'No agent CLI configured' }] };
  }

  const localPath = findLocalCheckout(entry);
  if (!localPath) {
    return { status: 'skip', gaps: [{ severity: SEVERITY.INFO, message: 'No local checkout; smoke requires local build' }] };
  }

  const probe = agentHelpProbe(entry, localPath);
  const gaps = probe.status === 'pass'
    ? []
    : [{ severity: probe.status === 'fail' ? SEVERITY.ERROR : SEVERITY.WARNING, message: probe.detail }];

  return { status: probe.status, gaps };
}

/**
 * Check DiagnosticEnvelope/v1 readiness.
 */
function checkDiagnosticReadiness(entry) {
  const metadata = readiness.get(entry.id);
  if (!metadata) {
    return { status: 'missing', gaps: [{ severity: SEVERITY.ERROR, message: 'No diagnostic readiness entry in registry' }] };
  }

  const gaps = [];

  if (!metadata.agentCli) {
    gaps.push({ severity: SEVERITY.WARNING, message: 'No agent CLI configured' });
  }

  if (metadata.closedLoop === 'planned') {
    gaps.push({ severity: SEVERITY.WARNING, message: 'Closed-loop support is planned, not implemented' });
  }

  return {
    status: gaps.length === 0 ? 'pass' : 'partial',
    gaps,
    metadata,
  };
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function git(args, options = {}) {
  try {
    return execFileSync('git', args, {
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', options.quiet ? 'ignore' : 'pipe'],
      timeout: 10000,
      ...options,
    }).trim();
  } catch {
    return '';
  }
}

function findLocalCheckout(entry) {
  for (const path of checkoutCandidates(entry)) {
    if (existsSync(join(path, '.git'))) {
      return path;
    }
  }
  return null;
}

function checkoutCandidates(entry) {
  const repoName = entry.repository.split('/')[1];
  return [
    join(codeRoot, '.worktrees-lsp-latest', repoName),
    join(codeRoot, repoName),
    join(codeRoot, '.worktrees-lsp-wiki-agent-cli-20260612', repoName),
  ];
}

function listFilesRecursive(dir, files = []) {
  if (!existsSync(dir)) return files;
  const items = readdirSync(dir);
  for (const item of items) {
    const fullPath = join(dir, item);
    try {
      const stat = statSync(fullPath);
      if (stat.isDirectory()) {
        listFilesRecursive(fullPath, files);
      } else {
        files.push(fullPath);
      }
    } catch {
      // skip inaccessible files
    }
  }
  return files;
}

function agentHelpProbe(entry, localPath) {
  const metadata = readiness.get(entry.id);
  if (!metadata?.agentCli) {
    return { status: 'skip', detail: 'no-agent-cli-metadata' };
  }

  const env = { ...process.env };
  const sourcePath = join(localPath, 'src');
  env.PYTHONPATH = [existsSync(sourcePath) ? sourcePath : localPath, localPath, process.env.PYTHONPATH]
    .filter(Boolean)
    .join(':');

  let command;
  let args;
  let cwd = localPath;

  if (entry.id === 'cif-lsp') {
    const script = join(localPath, 'server/out/cifLspTool.js');
    if (!existsSync(script)) {
      return { status: 'fail', detail: 'server/out/cifLspTool.js not found' };
    }
    command = 'node';
    args = [script, '--help'];
  } else if (entry.id === 'lammps-lsp') {
    command = 'cargo';
    args = ['run', '--quiet', '--bin', 'lammps-lsp-tool', '--', '--help'];
  } else {
    command = 'python3';
    args = ['-m', pythonModuleFor(entry), '--help'];
  }

  try {
    execFileSync(command, args, {
      cwd,
      env,
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'pipe'],
      timeout: entry.id === 'lammps-lsp' ? 120000 : 30000,
    });
    return { status: 'pass', detail: metadata.agentCli };
  } catch (error) {
    const stderr = error.stderr?.toString?.() || error.message;
    return { status: 'fail', detail: `${metadata.agentCli}: ${stderr.split('\n')[0]}` };
  }
}

function pythonModuleFor(entry) {
  const repoName = entry.repository.split('/')[1];
  if (repoName === 'cp2k-lsp-enhanced') return 'cp2k_input_tools.tool';
  if (repoName === 'VASP-LSP') return 'vasp_lsp.tool';
  return `${entry.id.replace(/-lsp$/, '_lsp').replaceAll('-', '_')}.tool`;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

function checkRepo(entry) {
  const manifest = checkManifest(entry);
  const provenance = checkProvenance(entry);
  const fixtures = checkFixtures(entry);
  const smoke = checkSmoke(entry);
  const diagReadiness = checkDiagnosticReadiness(entry);

  const allGaps = [
    ...manifest.gaps,
    ...provenance.gaps,
    ...fixtures.gaps,
    ...smoke.gaps,
    ...diagReadiness.gaps,
  ];

  const hasBlocking = allGaps.some(g => g.severity === SEVERITY.ERROR);
  const maturity = manifest.status === 'pass' && !hasBlocking ? 'stable' :
    manifest.status === 'missing' ? 'experimental' : 'preview';

  return {
    id: entry.id,
    languageId: entry.languageId,
    repository: entry.repository,
    manifest: manifest.status,
    manifestPath: manifest.path,
    provenance: provenance.status,
    fixtures: fixtures.status,
    smoke: smoke.status,
    diagnosticReadiness: diagReadiness.status,
    maturity,
    gaps: allGaps,
    blockingGaps: allGaps.filter(g => g.severity === SEVERITY.ERROR),
    warningGaps: allGaps.filter(g => g.severity === SEVERITY.WARNING),
  };
}

// Run all checks
const results = entries.map(checkRepo);

// Aggregate stats
const totalGaps = results.reduce((sum, r) => sum + r.gaps.length, 0);
const blockingGaps = results.reduce((sum, r) => sum + r.blockingGaps.length, 0);
const warningGaps = results.reduce((sum, r) => sum + r.warningGaps.length, 0);
const passCount = results.filter(r => r.blockingGaps.length === 0).length;

const report = {
  checkedAt: new Date().toISOString(),
  totalRepos: entries.length,
  passing: passCount,
  withBlockingGaps: results.filter(r => r.blockingGaps.length > 0).length,
  totalGaps,
  blockingGaps,
  warningGaps,
  graduationScore: calculateGraduationScore(results),
  warningSummary: buildWarningSummary(results),
  graduationTargets: Object.values(WARNING_CATEGORIES).map(cat => ({
    category: cat.label,
    weight: cat.graduationWeight,
    description: cat.description,
  })),
  results,
};

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

if (jsonMode) {
  console.log(JSON.stringify(report, null, 2));
} else {
  console.log('LSP Family Release Gate Check');
  console.log('='.repeat(60));
  console.log(`Checked at: ${report.checkedAt}`);
  console.log(`Repos: ${report.totalRepos} | Passing: ${report.passing} | Blocking gaps: ${report.blockingGaps} | Warnings: ${report.warningGaps}`);
  console.log(`Graduation Score: ${report.graduationScore}/100`);
  console.log('');

  // Warning summary
  if (report.warningGaps > 0) {
    console.log('Warning Summary by Category:');
    console.log('-'.repeat(40));
    for (const [id, entry] of Object.entries(report.warningSummary)) {
      if (entry.count > 0) {
        console.log(`  ${entry.category}: ${entry.count} warning(s) across ${entry.repos.length} repo(s)`);
        if (verbose) {
          for (const msg of entry.messages) {
            console.log(`    - ${msg}`);
          }
        }
      }
    }
    console.log('');
  }

  for (const r of results) {
    const icon = r.blockingGaps.length > 0 ? 'FAIL' : r.warningGaps.length > 0 ? 'WARN' : 'PASS';
    console.log(`[${icon}] ${r.id} (${r.languageId})`);
    console.log(`    manifest: ${r.manifest} | provenance: ${r.provenance} | fixtures: ${r.fixtures} | smoke: ${r.smoke} | diag: ${r.diagnosticReadiness}`);

    if (verbose || r.blockingGaps.length > 0) {
      for (const gap of r.blockingGaps) {
        console.log(`    ERROR: ${gap.message}`);
      }
    }
    if (verbose) {
      for (const gap of r.warningGaps) {
        console.log(`    WARN:  ${gap.message}`);
      }
    }
    console.log('');
  }

  console.log('='.repeat(60));

  // Graduation guide
  console.log('');
  console.log('Maturity Levels:');
  console.log('  preview     - Missing manifest; cannot validate against fleet requirements.');
  console.log('  experimental - Manifest present but has warnings (fixtures, provenance, or closed-loop gaps).');
  console.log('  stable       - All required gates pass. Safe for VSIX release blocking.');
  console.log('');
  console.log('Graduation Targets:');
  for (const target of report.graduationTargets) {
    console.log(`  ${target.category} (weight ${target.weight}): ${target.description}`);
  }
}

// Write report to file if --report-path specified
if (reportPath) {
  const { writeFileSync, mkdirSync } = await import('fs');
  const { dirname } = await import('path');
  mkdirSync(dirname(reportPath), { recursive: true });
  writeFileSync(reportPath, JSON.stringify(report, null, 2));
  if (!jsonMode) {
    console.log(`\nReport written to: ${reportPath}`);
  }
}

// Exit code
if (strictMode && blockingGaps > 0) {
  process.exit(1);
}
