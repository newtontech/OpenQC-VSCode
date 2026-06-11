#!/usr/bin/env node

import { execFileSync } from 'child_process';
import { existsSync } from 'fs';
import { dirname, join, resolve } from 'path';
import { fileURLToPath } from 'url';
import { readFileSync } from 'fs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '..');
const codeRoot = resolve(repoRoot, '..');
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
  const primary = readLocalCheckout(join(codeRoot, repoName));
  const latestWorktree = readLocalCheckout(join(codeRoot, '.worktrees-lsp-latest', repoName));

  if (primary.head === remote) {
    return { ...primary, source: 'sibling' };
  }

  if (latestWorktree.head === remote) {
    return { ...latestWorktree, source: 'latest-worktree' };
  }

  return { ...primary, source: 'sibling' };
}

const rows = entries.map(entry => {
  const remote = remoteHead(entry);
  const local = localHead(entry, remote);
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
    remote: remote.slice(0, 12),
    local: local.head.slice(0, 12) || '-',
    source: local.source,
    dirty: local.dirty ? 'yes' : 'no',
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
  status: row.status,
})));

const drift = rows.filter(row => row.status === 'not-at-remote-head' || row.status === 'missing-local');
if (drift.length > 0) {
  console.log(`\n${drift.length} LSP checkout(s) are not at the configured remote default-branch HEAD.`);
  console.log('Run again after updating sibling checkouts, or use this as release-note evidence.');
}

if (failOnDrift && drift.length > 0) {
  process.exitCode = 1;
}
