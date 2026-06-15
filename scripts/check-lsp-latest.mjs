#!/usr/bin/env node

import { execFileSync } from 'child_process';
import { existsSync } from 'fs';
import { basename, dirname, join, resolve } from 'path';
import { fileURLToPath } from 'url';
import { readFileSync } from 'fs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '..');
const codeRoot = resolveCodeRoot(repoRoot);
const registryPath = join(repoRoot, 'src/lsp/registry.ts');
const failOnDrift = process.argv.includes('--fail-on-drift');

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

function git(args, options = {}) {
  return execFileSync('git', args, {
    encoding: 'utf8',
    stdio: ['ignore', 'pipe', options.quiet ? 'ignore' : 'pipe'],
    ...options,
  }).trim();
}

function resolveCodeRoot(root) {
  const parent = resolve(root, '..');
  return basename(parent) === '.worktrees' ? resolve(parent, '..') : parent;
}

function remoteHead(entry) {
  const url = `https://github.com/${entry.repository}.git`;
  const output = git(['ls-remote', url, `refs/heads/${entry.defaultBranch}`]);
  return output.split(/\s+/)[0] || '';
}

function readLocalCheckout(localPath) {
  if (!existsSync(join(localPath, '.git'))) {
    return { path: localPath, head: '', dirty: '', exists: false };
  }

  const head = git(['-C', localPath, 'rev-parse', 'HEAD'], { quiet: true });
  const dirty = git(['-C', localPath, 'status', '--short'], { quiet: true });
  return { path: localPath, head, dirty, exists: true };
}

function localHead(entry, remote) {
  const repoName = entry.repository.split('/')[1];
  const candidates = [
    { source: 'sibling', checkout: readLocalCheckout(join(codeRoot, repoName)) },
    { source: 'latest-worktree', checkout: readLocalCheckout(join(codeRoot, '.worktrees-lsp-latest', repoName)) },
    {
      source: 'wiki-agent-worktree',
      checkout: readLocalCheckout(join(codeRoot, '.worktrees-lsp-wiki-agent-cli-20260612', repoName)),
    },
  ];

  const primary = candidates[0].checkout;

  for (const candidate of candidates) {
    if (candidate.checkout.head === remote) {
      return { ...candidate.checkout, source: candidate.source };
    }
  }

  return { ...primary, source: 'sibling' };
}

function pythonModuleFor(entry) {
  const repoName = entry.repository.split('/')[1];
  if (repoName === 'cp2k-lsp-enhanced') return 'cp2k_input_tools.tool';
  if (repoName === 'VASP-LSP') return 'vasp_lsp.tool';
  const packageName = entry.id.replace(/-lsp$/, '_lsp').replaceAll('-', '_');
  return `${packageName}.tool`;
}

function agentHelpProbe(entry, local) {
  const metadata = readiness.get(entry.id);
  if (!metadata?.agentCli) {
    return { status: 'unavailable', detail: 'no-agent-cli-metadata' };
  }
  if (!local.exists) {
    return { status: 'skipped', detail: 'missing-local-checkout' };
  }

  const repoName = entry.repository.split('/')[1];
  const localPath = local.path;
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
      return { status: 'missing-build', detail: 'server/out/cifLspTool.js not found' };
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

const rows = entries.map(entry => {
  const remote = remoteHead(entry);
  const local = localHead(entry, remote);
  const metadata = readiness.get(entry.id);
  const agentHelp = agentHelpProbe(entry, local);
  const status = !local.exists
    ? 'missing-local'
    : local.head === remote
      ? local.dirty
        ? 'latest-with-local-changes'
        : local.source === 'latest-worktree'
          ? 'latest-via-worktree'
          : 'latest'
      : 'not-at-remote-head';

  return {
    ...entry,
    agentCli: metadata?.agentCli ?? '-',
    remote: remote.slice(0, 12),
    local: local.head.slice(0, 12) || '-',
    source: local.source,
    dirty: local.dirty ? 'yes' : 'no',
    agentHelp: agentHelp.status,
    agentHelpDetail: agentHelp.detail,
    status,
  };
});

console.table(rows.map(row => ({
  id: row.id,
  language: row.languageId,
  branch: row.defaultBranch,
  source: row.source,
  local: row.local,
  remote: row.remote,
  dirty: row.dirty,
  agentCli: row.agentCli,
  agentHelp: row.agentHelp,
  status: row.status,
})));

const drift = rows.filter(row => row.status === 'not-at-remote-head' || row.status === 'missing-local');
const agentFailures = rows.filter(row => row.agentHelp === 'fail' || row.agentHelp === 'missing-build');
if (drift.length > 0) {
  console.log(`\n${drift.length} LSP checkout(s) are not at the configured remote default-branch HEAD.`);
  console.log('Run again after updating sibling checkouts, or use this as release-note evidence.');
}

if (agentFailures.length > 0) {
  console.log(`\n${agentFailures.length} LSP agent CLI probe(s) failed.`);
  for (const failure of agentFailures) {
    console.log(`- ${failure.id}: ${failure.agentHelpDetail}`);
  }
}

if (failOnDrift && (drift.length > 0 || agentFailures.length > 0)) {
  process.exitCode = 1;
}
