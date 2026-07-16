#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const scriptRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const root = process.env.OPENQC_RELEASE_ROOT
  ? path.resolve(process.env.OPENQC_RELEASE_ROOT)
  : scriptRoot;
const packageJson = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
const expectedTag = `v${packageJson.version}`;
const args = process.argv.slice(2);
const prepare = args.includes('--prepare');
const noFetch = args.includes('--no-fetch');
const tagIndex = args.indexOf('--tag');
const tag = tagIndex >= 0 ? args[tagIndex + 1] : expectedTag;

function git(gitArgs, { allowFailure = false } = {}) {
  const result = spawnSync('git', gitArgs, { cwd: root, encoding: 'utf8' });
  if (result.status !== 0 && !allowFailure) {
    process.stderr.write(result.stderr || result.stdout || 'git command failed\n');
    process.exit(result.status ?? 1);
  }
  return { status: result.status ?? 1, stdout: result.stdout.trim() };
}

function fail(message) {
  console.error(`Release guard failed: ${message}`);
  process.exit(1);
}

if (!tag || tag !== expectedTag) {
  fail(`tag ${tag || '<missing>'} must match package.json version as ${expectedTag}`);
}

if (!noFetch) {
  git(['fetch', '--quiet', 'origin', 'master', '--tags']);
}

const status = git(['status', '--porcelain', '--untracked-files=normal']).stdout;
if (status) {
  fail(`working tree is not clean:\n${status}`);
}

const head = git(['rev-parse', 'HEAD']).stdout;
const originMaster = git(['rev-parse', 'refs/remotes/origin/master']).stdout;
if (head !== originMaster) {
  fail(`HEAD ${head} is not synchronized with origin/master ${originMaster}`);
}

const tagRef = `refs/tags/${tag}`;
const tagExists =
  git(['show-ref', '--verify', '--quiet', tagRef], { allowFailure: true }).status === 0;

if (prepare) {
  if (tagExists) {
    fail(`tag already exists: ${tag}`);
  }
  console.log(`Release tag ${tag} may be created at synchronized commit ${head}.`);
  process.exit(0);
}

if (!tagExists) {
  fail(`tag does not exist in the checkout: ${tag}`);
}

const tagCommit = git(['rev-parse', `${tag}^{commit}`]).stdout;
if (tagCommit !== head) {
  fail(`tag ${tag} points to ${tagCommit}, not HEAD ${head}`);
}

console.log(`Verified ${tag} at clean origin/master commit ${head}.`);
