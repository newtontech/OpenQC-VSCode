#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const packageJson = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
const defaultVsix = path.join(root, 'vsix', `${packageJson.name}-${packageJson.version}.vsix`);
const vsixPath = path.resolve(process.argv[2] || defaultVsix);

function unzip(args) {
  const result = spawnSync('unzip', args, { cwd: root, encoding: 'utf8' });
  if (result.status !== 0) {
    process.stderr.write(result.stderr || result.stdout || 'unzip failed\n');
    process.exit(result.status ?? 1);
  }
  return result.stdout;
}

if (!fs.existsSync(vsixPath)) {
  throw new Error(`VSIX not found: ${vsixPath}. Run \`npm run package:vsix\` first.`);
}

const entries = unzip(['-Z1', vsixPath])
  .split(/\r?\n/)
  .map(entry => entry.replace(/\\/g, '/'))
  .filter(Boolean);

const forbidden = [
  /^extension\/(?:\.git|\.github|\.vscode|\.vscode-test|\.worktrees|\.omo)(?:\/|$)/,
  /^extension\/(?:ai-protocols|core|coverage|docs|output|screenshots|scripts|server|tests|vsix)(?:\/|$)/,
  /^extension\/(?:AGENTS|MEMORY)\.md$/,
  /(?:^|\/)\.env(?:\.|$)/,
  /(?:^|\/)(?:id_rsa|credentials\.json)$/,
  /\.(?:key|pem|pyc|vsix|DS_Store)$/,
  /(?:^|\/)__pycache__(?:\/|$)/,
];
const forbiddenEntries = entries.filter(entry => forbidden.some(pattern => pattern.test(entry)));

if (forbiddenEntries.length > 0) {
  throw new Error(`VSIX contains forbidden files:\n${forbiddenEntries.join('\n')}`);
}

const requiredEntries = [
  'extension/out/extension.js',
  'extension/python/format_converter.py',
  'extension/node_modules/3dmol/build/3Dmol-min.js',
  'extension/node_modules/ngl/dist/ngl.js',
  'extension/node_modules/plotly.js-dist-min/plotly.min.js',
  'extension/node_modules/vscode-languageclient/lib/node/main.js',
];
const missingEntries = requiredEntries.filter(entry => !entries.includes(entry));
if (missingEntries.length > 0) {
  throw new Error(`VSIX is missing required runtime files:\n${missingEntries.join('\n')}`);
}

const manifest = JSON.parse(unzip(['-p', vsixPath, 'extension/package.json']));
const expectedIdentity = `${packageJson.publisher}.${packageJson.name}@${packageJson.version}`;
const actualIdentity = `${manifest.publisher}.${manifest.name}@${manifest.version}`;

if (actualIdentity !== expectedIdentity) {
  throw new Error(`VSIX identity mismatch: expected ${expectedIdentity}, got ${actualIdentity}`);
}

if (manifest.displayName !== packageJson.displayName) {
  throw new Error(
    `VSIX display name mismatch: expected ${packageJson.displayName}, got ${manifest.displayName}`
  );
}

const stat = fs.statSync(vsixPath);
const maxBytes = 25 * 1024 * 1024;
const maxEntries = 900;
if (stat.size > maxBytes) {
  throw new Error(`VSIX is unexpectedly large: ${stat.size} bytes (limit ${maxBytes})`);
}
if (entries.length > maxEntries) {
  throw new Error(`VSIX has too many files: ${entries.length} (limit ${maxEntries})`);
}

console.log(
  `Verified ${expectedIdentity}: ${entries.length} files, ${(stat.size / 1024 / 1024).toFixed(2)} MiB`
);
