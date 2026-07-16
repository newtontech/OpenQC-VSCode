#!/usr/bin/env node

import { spawnSync } from 'child_process';
import fs from 'fs';
import path from 'path';

const root = path.resolve(new URL('..', import.meta.url).pathname);
const packageJson = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
const defaultVsix = path.join(root, 'vsix', `${packageJson.name}-${packageJson.version}.vsix`);
const vsixPath = path.resolve(process.argv[2] || defaultVsix);
const extractRoot = path.join(root, 'output', 'playwright', 'vsix-visual-smoke');
const extensionRoot = path.join(extractRoot, 'extension');

if (!fs.existsSync(vsixPath)) {
  throw new Error(`VSIX not found: ${vsixPath}. Run \`npm run package:vsix\` first.`);
}

fs.rmSync(extractRoot, { recursive: true, force: true });
fs.mkdirSync(extractRoot, { recursive: true });

const unzip = spawnSync('unzip', ['-q', vsixPath, '-d', extractRoot], {
  cwd: root,
  encoding: 'utf8',
});

if (unzip.status !== 0) {
  process.stderr.write(unzip.stdout || '');
  process.stderr.write(unzip.stderr || '');
  throw new Error(`Failed to extract VSIX for visual smoke: ${vsixPath}`);
}

const requiredFiles = [
  path.join(extensionRoot, 'out', 'webviews', 'openqcViewerWebview.js'),
  path.join(extensionRoot, 'media', 'openqc-viewer.js'),
  path.join(extensionRoot, 'media', 'openqc-viewer.css'),
  path.join(extensionRoot, 'media', 'vendor', '3dmol', '3Dmol-min.js'),
];

for (const file of requiredFiles) {
  if (!fs.existsSync(file)) {
    throw new Error(`Extracted VSIX is missing required visual runtime file: ${file}`);
  }
}

const smoke = spawnSync(
  process.execPath,
  [
    path.join(root, 'scripts', 'visual-smoke-openqc-viewer.mjs'),
    '--extension-root',
    extensionRoot,
    '--output-prefix',
    'openqc-viewer-vsix-smoke',
  ],
  { cwd: root, encoding: 'utf8', stdio: 'inherit' }
);

if (smoke.status !== 0) {
  throw new Error(`Packaged VSIX visual smoke failed for ${vsixPath}`);
}
