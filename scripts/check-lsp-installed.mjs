#!/usr/bin/env node

import { execFileSync } from 'child_process';
import { existsSync, readFileSync, writeFileSync } from 'fs';
import { basename, dirname, join, resolve } from 'path';
import { fileURLToPath } from 'url';

const scriptDir = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(scriptDir, '..');
const codeRoot = process.env.OPENQC_LSP_CODE_ROOT
  ? resolve(process.env.OPENQC_LSP_CODE_ROOT)
  : resolveCodeRoot(repoRoot);
const runtimeRoot = process.env.OPENQC_LSP_RUNTIME_ROOT
  ? resolve(process.env.OPENQC_LSP_RUNTIME_ROOT)
  : join(codeRoot, '.lsp-latest');
const registryPath = process.env.OPENQC_LSP_REGISTRY_PATH
  ? resolve(process.env.OPENQC_LSP_REGISTRY_PATH)
  : join(repoRoot, 'src/lsp/registry.ts');
const matrixPath = process.env.OPENQC_LSP_RELEASE_MATRIX_PATH
  ? resolve(process.env.OPENQC_LSP_RELEASE_MATRIX_PATH)
  : join(scriptDir, 'lsp-release-matrix.json');
const strict = process.argv.includes('--strict') || process.argv.includes('--fail-on-drift');
const jsonMode = process.argv.includes('--json');
const reportPath = process.argv.find((arg, index) => process.argv[index - 1] === '--report-path');

function resolveCodeRoot(root) {
  const parent = resolve(root, '..');
  if (basename(parent) !== '.worktrees') return parent;
  const worktreeContainer = resolve(parent, '..');
  return existsSync(join(worktreeContainer, '.git'))
    ? resolve(worktreeContainer, '..')
    : worktreeContainer;
}

function gitHead(path) {
  if (!existsSync(path)) return '';
  try {
    return execFileSync('git', ['-C', path, 'rev-parse', 'HEAD'], {
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'ignore'],
    }).trim();
  } catch {
    return '';
  }
}

function headTag(path) {
  if (!existsSync(path)) return '';
  try {
    return (
      execFileSync('git', ['-C', path, 'tag', '--points-at', 'HEAD', '--sort=-version:refname'], {
        encoding: 'utf8',
        stdio: ['ignore', 'pipe', 'ignore'],
      })
        .trim()
        .split('\n')
        .filter(Boolean)[0] ?? ''
    );
  } catch {
    return '';
  }
}

function manifestVersion(path) {
  const packageJson = join(path, 'package.json');
  if (existsSync(packageJson)) {
    try {
      return JSON.parse(readFileSync(packageJson, 'utf8')).version ?? '';
    } catch {
      return '';
    }
  }

  for (const name of ['pyproject.toml', 'Cargo.toml']) {
    const manifest = join(path, name);
    if (!existsSync(manifest)) continue;
    const content = readFileSync(manifest, 'utf8');
    const projectSection =
      content.match(/\[(?:project|package)\]([\s\S]*?)(?=\n\[|$)/)?.[1] ?? content;
    const match = projectSection.match(/^version\s*=\s*["']([^"']+)["']/m);
    if (match) return match[1];
  }

  const versionFile = join(path, 'VERSION');
  return existsSync(versionFile) ? readFileSync(versionFile, 'utf8').trim() : '';
}

function sourceCheckout(repoName) {
  const candidates = [join(codeRoot, '.worktrees-lsp-latest', repoName), join(codeRoot, repoName)];
  return candidates.find(candidate => gitHead(candidate)) ?? candidates[0];
}

function executablePaths(entry, agentCli) {
  const repoName = entry.repository.split('/')[1];
  const installedCheckout = join(runtimeRoot, repoName);
  if (entry.id === 'cif-lsp') {
    return {
      server: join(installedCheckout, 'server', 'out', 'server.js'),
      agent: join(installedCheckout, 'server', 'out', 'cifLspTool.js'),
    };
  }
  if (entry.id === 'lammps-lsp') {
    return {
      server: join(installedCheckout, 'target', 'release', entry.executable),
      agent: join(installedCheckout, 'target', 'release', agentCli),
    };
  }
  return {
    server: join(runtimeRoot, '.venv', 'bin', entry.executable),
    agent: join(runtimeRoot, '.venv', 'bin', agentCli),
  };
}

const registry = readFileSync(registryPath, 'utf8');
const entries = [
  ...registry.matchAll(
    /id: '([^']+)'[\s\S]*?repository: '([^']+)'[\s\S]*?executable: '([^']+)'[\s\S]*?languageId: '([^']+)'[\s\S]*?defaultBranch: '([^']+)'/g
  ),
].map(match => ({
  id: match[1],
  repository: match[2],
  executable: match[3],
  languageId: match[4],
  defaultBranch: match[5],
}));
const readiness = new Map(
  [...registry.matchAll(/'([^']+)': diagnosticReadiness\('([^']+)',/g)].map(match => [
    match[1],
    match[2],
  ])
);
const matrix = JSON.parse(readFileSync(matrixPath, 'utf8'));

if (entries.length !== 17) {
  throw new Error(`Expected 17 bundled LSP entries, found ${entries.length}`);
}
if (matrix.schemaVersion !== 'openqc.lsp.release-matrix.v1') {
  throw new Error(`Unsupported release matrix schema: ${matrix.schemaVersion}`);
}

const rows = entries.map(entry => {
  const repoName = entry.repository.split('/')[1];
  const sourcePath = sourceCheckout(repoName);
  const installedPath = join(runtimeRoot, repoName);
  const sourceCommit = gitHead(sourcePath);
  const installedCommit = gitHead(installedPath);
  const tagVersion = headTag(sourcePath).replace(/^v/, '');
  const sourceVersion = manifestVersion(sourcePath);
  const installedVersion = manifestVersion(installedPath);
  const release = matrix.entries[entry.id];
  const agentCli = readiness.get(entry.id) ?? '';
  const executables = executablePaths(entry, agentCli);
  const problems = [];

  if (!release) problems.push('missing-release-matrix-entry');
  if (!sourceCommit) problems.push('missing-source-checkout');
  if (!installedCommit) problems.push('missing-installed-checkout');
  if (sourceCommit && installedCommit && sourceCommit !== installedCommit) {
    problems.push('installed-commit-drift');
  }
  if (!existsSync(executables.server)) problems.push('missing-server-executable');
  if (!agentCli || !existsSync(executables.agent)) problems.push('missing-agent-cli');
  if (release?.targetVersion && installedVersion && release.targetVersion !== installedVersion) {
    problems.push('installed-version-drift');
  }

  return {
    id: entry.id,
    repository: entry.repository,
    languageId: entry.languageId,
    defaultBranch: entry.defaultBranch,
    package: release?.package ?? '',
    channel: release?.channel ?? '',
    sourcePath,
    installedPath,
    sourceCommit,
    installedCommit,
    tagVersion,
    manifestVersion: sourceVersion,
    installedVersion,
    registryVersion: release?.targetVersion ?? '',
    serverExecutable: executables.server,
    agentExecutable: executables.agent,
    runtimeParity: problems.length === 0 ? 'pass' : 'fail',
    problems,
  };
});

const report = {
  schemaVersion: 'openqc.lsp.runtime-ledger.v1',
  generatedAt: new Date().toISOString(),
  codeRoot,
  runtimeRoot,
  summary: {
    total: rows.length,
    passing: rows.filter(row => row.runtimeParity === 'pass').length,
    failing: rows.filter(row => row.runtimeParity === 'fail').length,
  },
  entries: rows,
};

if (reportPath) {
  writeFileSync(resolve(reportPath), `${JSON.stringify(report, null, 2)}\n`);
}

if (jsonMode) {
  process.stdout.write(`${JSON.stringify(report, null, 2)}\n`);
} else {
  console.table(
    rows.map(row => ({
      id: row.id,
      source: row.sourceCommit.slice(0, 12) || '-',
      installed: row.installedCommit.slice(0, 12) || '-',
      manifest: row.manifestVersion || '-',
      installedVersion: row.installedVersion || '-',
      target: row.registryVersion || '-',
      parity: row.runtimeParity,
      problems: row.problems.join(',') || '-',
    }))
  );
  console.log(
    `\nInstalled LSP runtime parity: ${report.summary.passing}/${report.summary.total} passing, ${report.summary.failing} failing.`
  );
}

if (strict && report.summary.failing > 0) {
  process.exitCode = 1;
}
