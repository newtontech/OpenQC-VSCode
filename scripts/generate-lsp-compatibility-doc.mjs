#!/usr/bin/env node

/**
 * Generate the OpenQC-facing LSP compatibility and maturity document from the
 * registry plus the family gate report.
 *
 * Usage:
 *   node scripts/generate-lsp-compatibility-doc.mjs
 *   node scripts/generate-lsp-compatibility-doc.mjs --family-report /tmp/report.json
 *   node scripts/generate-lsp-compatibility-doc.mjs --check
 */

import { execFileSync } from 'child_process';
import { existsSync, mkdirSync, readFileSync, writeFileSync } from 'fs';
import { basename, dirname, join, relative, resolve } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '..');
const codeRoot = resolveCodeRoot(repoRoot);

function resolveCodeRoot(root) {
  const parent = resolve(root, '..');
  if (basename(parent) !== '.worktrees') return parent;
  const worktreeContainer = resolve(parent, '..');
  return existsSync(join(worktreeContainer, '.git'))
    ? resolve(worktreeContainer, '..')
    : worktreeContainer;
}

const args = process.argv.slice(2);
const checkMode = args.includes('--check');
const dryRun = args.includes('--dry-run');
const outputPath = resolve(readFlag('--output') ?? join(repoRoot, 'docs', 'LSP_COMPATIBILITY.md'));
const familyReportPath = readFlag('--family-report');

const registryPath = join(repoRoot, 'src', 'lsp', 'registry.ts');
const registrySource = readFileSync(registryPath, 'utf8');
const registryEntries = parseRegistry(registrySource);
const readiness = parseReadiness(registrySource);
const familyReport = loadFamilyReport(familyReportPath);
const familyById = new Map((familyReport.results ?? []).map(result => [result.id, result]));

if (registryEntries.length === 0) {
  throw new Error(`No bundled LSP registry entries found in ${registryPath}`);
}

const markdown = renderDocument(registryEntries, readiness, familyById, familyReport);

if (checkMode) {
  const existing = existsSync(outputPath) ? readFileSync(outputPath, 'utf8') : '';
  if (existing !== markdown) {
    console.error(`${relative(repoRoot, outputPath)} is out of date. Run npm run lsp:generate-compatibility-doc.`);
    process.exit(1);
  }
  console.log(`${relative(repoRoot, outputPath)} is up to date.`);
} else if (dryRun) {
  process.stdout.write(markdown);
} else {
  mkdirSync(dirname(outputPath), { recursive: true });
  writeFileSync(outputPath, markdown);
  console.log(`Wrote ${relative(repoRoot, outputPath)} for ${registryEntries.length} LSP backends.`);
}

function readFlag(name) {
  const index = args.indexOf(name);
  if (index === -1) return null;
  const value = args[index + 1];
  if (!value || value.startsWith('--')) {
    throw new Error(`${name} requires a value`);
  }
  return value;
}

function loadFamilyReport(path) {
  if (path) {
    return JSON.parse(readFileSync(resolve(path), 'utf8'));
  }

  const output = execFileSync(process.execPath, [join(__dirname, 'check-lsp-family.mjs'), '--json'], {
    cwd: repoRoot,
    encoding: 'utf8',
    stdio: ['ignore', 'pipe', 'inherit'],
  });
  return JSON.parse(output);
}

function parseRegistry(source) {
  const marker = 'const BUNDLED_LSP_SERVERS';
  const start = source.indexOf(marker);
  if (start === -1) return [];

  const arrayStart = source.indexOf('[', start);
  const arrayEnd = source.indexOf('] as const', arrayStart);
  const arraySource = source.slice(arrayStart, arrayEnd);

  const blocks = [];
  let depth = 0;
  let blockStart = -1;
  for (let i = 0; i < arraySource.length; i++) {
    const char = arraySource[i];
    if (char === '{') {
      if (depth === 0) blockStart = i;
      depth++;
    } else if (char === '}') {
      depth--;
      if (depth === 0 && blockStart !== -1) {
        blocks.push(arraySource.slice(blockStart, i + 1));
        blockStart = -1;
      }
    }
  }

  return blocks.map(block => ({
    id: readStringField(block, 'id'),
    name: readStringField(block, 'name'),
    repository: readStringField(block, 'repository'),
    executable: readStringField(block, 'executable'),
    languageId: readStringField(block, 'languageId'),
    fileExtensions: readStringArray(block, 'fileExtensions'),
    fileNames: readStringArray(block, 'fileNames'),
    stability: readStringField(block, 'stability'),
    defaultBranch: readStringField(block, 'defaultBranch'),
  })).filter(entry => entry.id);
}

function parseReadiness(source) {
  const map = new Map();
  const pattern = /'([^']+)': diagnosticReadiness\('([^']+)', '([^']+)'(?:, '([^']+)')?\)/g;
  for (const match of source.matchAll(pattern)) {
    map.set(match[1], {
      agentCli: match[2],
      closedLoop: match[3],
      agentCliSmokeStatus: match[4] ?? 'pending',
    });
  }
  return map;
}

function readStringField(block, field) {
  const match = block.match(new RegExp(`\\b${field}:\\s*'([^']*)'`));
  return match?.[1] ?? '';
}

function readStringArray(block, field) {
  const start = block.search(new RegExp(`\\b${field}:\\s*\\[`));
  if (start === -1) return [];
  const open = block.indexOf('[', start);
  let depth = 0;
  for (let i = open; i < block.length; i++) {
    if (block[i] === '[') depth++;
    if (block[i] === ']') {
      depth--;
      if (depth === 0) {
        return [...block.slice(open + 1, i).matchAll(/'([^']+)'/g)].map(match => match[1]);
      }
    }
  }
  return [];
}

function renderDocument(entries, readiness, familyById, report) {
  const lines = [
    '# LSP Compatibility Matrix',
    '',
    'This document is generated from `src/lsp/registry.ts` and the LSP family gate report.',
    'Regenerate it with `npm run lsp:generate-compatibility-doc` after registry, backend, or gate changes.',
    '',
    '## Release Gate Summary',
    '',
    `| Metric | Value |`,
    `|--------|-------|`,
    `| Backends | ${report.totalRepos ?? entries.length} |`,
    `| Passing | ${report.passing ?? 'unknown'} |`,
    `| Blocking gaps | ${report.blockingGaps ?? 'unknown'} |`,
    `| Warnings | ${report.warningGaps ?? 'unknown'} |`,
    `| Graduation score | ${report.graduationScore ?? 'unknown'}/100 |`,
    '',
    '## Backend Matrix',
    '',
    '| Backend | Language | Branch | File Types | Registry Stability | Gate Maturity | Release Evidence | Manifest | Provenance | Fixtures | Smoke | Diagnostics | Traceability |',
    '|---------|----------|--------|------------|--------------------|---------------|------------------|----------|------------|----------|-------|-------------|--------------|',
  ];

  for (const entry of entries) {
    const family = familyById.get(entry.id);
    const evidence = family?.provenanceEvidence ?? buildLocalProvenanceEvidence(entry);
    lines.push([
      repoLink(entry),
      code(entry.languageId),
      code(entry.defaultBranch),
      formatFileTypes(entry),
      entry.stability,
      family?.maturity ?? 'unknown',
      formatReleaseEvidence(evidence),
      statusCell(family?.manifest),
      statusCell(family?.provenance),
      statusCell(family?.fixtures),
      statusCell(family?.smoke),
      statusCell(family?.diagnosticReadiness),
      statusCell(family?.traceability),
    ].join(' | ').replace(/^/, '| ').replace(/$/, ' |'));
  }

  lines.push(
    '',
    '## Agent CLI Readiness',
    '',
    '| Backend | Agent CLI | Operations | Help Smoke | Closed Loop |',
    '|---------|-----------|------------|------------|-------------|'
  );

  for (const entry of entries) {
    const ready = readiness.get(entry.id) ?? {};
    const family = familyById.get(entry.id);
    lines.push([
      code(entry.id),
      code(ready.agentCli ?? 'unconfigured'),
      '`check`, `context`, `complete`, `hover`, `symbols`, `fix`',
      family?.smoke ?? ready.agentCliSmokeStatus ?? 'unknown',
      ready.closedLoop ?? 'unknown',
    ].join(' | ').replace(/^/, '| ').replace(/$/, ' |'));
  }

  lines.push('', '## Actionable Gate Gaps', '');
  const gapLines = [];
  for (const entry of entries) {
    const family = familyById.get(entry.id);
    const gaps = [
      ...(family?.blockingGaps ?? []).map(gap => `ERROR: ${gap.message}`),
      ...(family?.warningGaps ?? []).map(gap => `WARN: ${gap.message}`),
    ];
    const actions = family?.provenanceActions ?? [];
    if (gaps.length === 0 && actions.length === 0) continue;
    gapLines.push(`### ${entry.id}`);
    for (const gap of gaps) gapLines.push(`- ${gap}`);
    for (const action of actions) gapLines.push(`- Action: ${action}`);
    gapLines.push('');
  }
  if (gapLines.length === 0) {
    lines.push('No blocking or warning gaps are reported by `scripts/check-lsp-family.mjs`.');
  } else {
    lines.push(...gapLines);
  }

  lines.push(
    '',
    '## OpenQC Integration Guarantees',
    '',
    '| Capability | OpenQC guarantee |',
    '|------------|------------------|',
    '| Language contribution | Every registry entry has a matching `contributes.languages` item in `package.json`. |',
    '| Grammar contribution | Every language has a TextMate grammar under `syntaxes/`. |',
    '| Configuration | Every language exposes `enabled`, `path`, `command`, `args`, and `env` settings under `openqc.lsp.<languageId>.*`. |',
    '| Startup | OpenQC starts the configured command over stdio and can prefer sibling local repositories when no user override is set. |',
    '| Latest tracking | `npm run lsp:check-latest -- --fail-on-drift` verifies local checkout cleanliness, remote HEAD parity, and agent CLI help. |',
    '| Family maturity | `npm run lsp:check-family -- --strict` verifies manifest, provenance, fixture, smoke, and Diagnostic Engine readiness. |',
    '| Traceability aggregation | `npm run lsp:check-family -- --strict` consumes backend docstring/wiki/raw traceability reports when present and reports missing or broken evidence links. |',
    '| Bohrium routing | `npm run lsp:check-bohrium-registry` verifies OpenQC and Bohrium backend ids and agent CLI names stay aligned. |',
    '',
    '## Known Limits',
    '',
    '- This gate proves OpenQC integration readiness and the shared agent-facing contract; it does not prove exhaustive grammar coverage for every historical software release.',
    '- Complete version-aware parser/rule coverage and official-doc ingestion remain tracked in each backend repository issue queue.',
    '- Runtime scientific correctness still requires backend-specific fixture suites and, where available, real executable/log smoke tests.',
    '- OpenQC launches the configured executable or sibling checkout; it does not vendor or pin standalone backend binaries.'
  );

  return `${lines.join('\n')}\n`;
}

function repoLink(entry) {
  return `[${entry.repository}](https://github.com/${entry.repository})`;
}

function code(value) {
  return `\`${value}\``;
}

function formatFileTypes(entry) {
  const extensions = entry.fileExtensions.map(ext => code(ext.startsWith('.') ? ext : `.${ext}`));
  const names = entry.fileNames.map(name => code(name));
  const values = [...extensions, ...names];
  return values.length > 0 ? values.join(', ') : 'none';
}

function statusCell(status) {
  if (!status) return 'unknown';
  return status === 'pass' ? 'pass' : status;
}

function formatReleaseEvidence(evidence) {
  const tag = evidence?.latestTag;
  const versionFile = firstExistingEvidencePath(evidence?.versionFiles);
  const version = versionFile ? readFirstLine(versionFile) : null;
  const manifestVersion = evidence?.manifestReleaseVersion;
  const head = evidence?.head ? String(evidence.head).slice(0, 12) : null;
  const parts = [];
  if (tag) {
    const remoteState = evidence?.remoteTagPresent === false ? ' local-only' : '';
    parts.push(`${code(tag)}${remoteState}`);
  }
  if (version) parts.push(`VERSION ${code(version)}`);
  if (!tag && manifestVersion) parts.push(`manifest ${code(manifestVersion)}`);
  if (head) parts.push(`HEAD ${code(head)}`);
  return parts.length > 0 ? parts.join('<br>') : 'unreleased/local';
}

function firstExistingEvidencePath(paths = []) {
  return paths.find(path => typeof path === 'string' && existsSync(path)) ?? null;
}

function readFirstLine(path) {
  try {
    return readFileSync(path, 'utf8').split(/\r?\n/).find(Boolean)?.trim() ?? null;
  } catch {
    return null;
  }
}

function buildLocalProvenanceEvidence(entry) {
  const localPath = findLocalCheckout(entry);
  if (!localPath) {
    return {};
  }
  const latestTag = git(['-C', localPath, 'tag', '--sort=-version:refname']).split(/\r?\n/).find(Boolean) ?? null;
  return {
    localPath,
    head: git(['-C', localPath, 'rev-parse', 'HEAD']),
    latestTag,
    remoteTagPresent: hasRemoteTag(localPath, latestTag),
    versionFiles: ['VERSION', 'version.txt'].map(file => join(localPath, file)).filter(path => existsSync(path)),
    changelogPaths: ['CHANGELOG.md', 'CHANGELOG'].map(file => join(localPath, file)).filter(path => existsSync(path)),
    manifestPath: existsSync(join(localPath, 'lsp-capabilities.json')) ? join(localPath, 'lsp-capabilities.json') : null,
  };
}

function findLocalCheckout(entry) {
  const repoName = entry.repository.split('/')[1];
  const candidates = [
    join(codeRoot, '.worktrees-lsp-latest', repoName),
    join(codeRoot, repoName),
    join(codeRoot, '.worktrees-lsp-wiki-agent-cli-20260612', repoName),
  ];
  return candidates.find(path => existsSync(join(path, '.git'))) ?? null;
}

function git(args) {
  try {
    return execFileSync('git', args, {
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'ignore'],
      timeout: 10000,
    }).trim();
  } catch {
    return '';
  }
}

function hasRemoteTag(localPath, tag) {
  if (!tag) return null;
  return git(['-C', localPath, 'ls-remote', '--tags', 'origin', `refs/tags/${tag}`]).length > 0;
}
