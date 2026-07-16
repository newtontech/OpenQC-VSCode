#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const packageJson = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
const tag = `v${packageJson.version}`;
const shouldPush = process.argv.includes('--push');

function run(command, args) {
  const result = spawnSync(command, args, { cwd: root, stdio: 'inherit' });
  if (result.status !== 0) {
    process.exit(result.status ?? 1);
  }
}

run(process.execPath, [path.join(root, 'scripts', 'release-guard.mjs'), '--prepare']);
run('git', ['tag', '-a', tag, '-m', `Release ${tag}`]);

if (shouldPush) {
  run('git', ['push', 'origin', tag]);
}

console.log(`Created release tag ${tag}${shouldPush ? ' and pushed it' : ''}.`);
