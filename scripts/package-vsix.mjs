#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const packageJson = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
const outDir = path.join(root, 'vsix');
const outPath = path.join(outDir, `${packageJson.name}-${packageJson.version}.vsix`);
const vsce = path.join(
  root,
  'node_modules',
  '.bin',
  process.platform === 'win32' ? 'vsce.cmd' : 'vsce'
);

if (!fs.existsSync(vsce)) {
  throw new Error('Pinned @vscode/vsce is not installed. Run `npm ci` first.');
}

fs.mkdirSync(outDir, { recursive: true });
fs.rmSync(outPath, { force: true });

const result = spawnSync(vsce, ['package', '--out', outPath], {
  cwd: root,
  stdio: 'inherit',
});

if (result.status !== 0) {
  process.exit(result.status ?? 1);
}

console.log(`Packaged VSIX: ${path.relative(root, outPath)}`);
