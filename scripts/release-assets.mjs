#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import crypto from 'node:crypto';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const packageJson = JSON.parse(fs.readFileSync(path.join(root, 'package.json'), 'utf8'));
const outDir = path.join(root, 'vsix');
const baseName = `${packageJson.name}-${packageJson.version}`;
const vsixPath = path.join(outDir, `${baseName}.vsix`);
const checksumPath = `${vsixPath}.sha256`;
const sbomPath = path.join(outDir, `${baseName}.sbom.json`);

if (!fs.existsSync(vsixPath)) {
  throw new Error(`VSIX not found: ${vsixPath}. Run \`npm run package:vsix\` first.`);
}

const digest = crypto.createHash('sha256').update(fs.readFileSync(vsixPath)).digest('hex');
fs.writeFileSync(checksumPath, `${digest}  ${path.basename(vsixPath)}\n`, 'utf8');

const sbom = spawnSync('npm', ['sbom', '--omit=dev', '--sbom-format=cyclonedx'], {
  cwd: root,
  encoding: 'utf8',
});
if (sbom.status !== 0) {
  process.stderr.write(sbom.stderr || sbom.stdout || 'npm sbom failed\n');
  process.exit(sbom.status ?? 1);
}

const parsedSbom = JSON.parse(sbom.stdout);
if (parsedSbom.bomFormat !== 'CycloneDX') {
  throw new Error(`Unexpected SBOM format: ${parsedSbom.bomFormat || 'missing'}`);
}
fs.writeFileSync(sbomPath, `${JSON.stringify(parsedSbom, null, 2)}\n`, 'utf8');

console.log(`SHA-256: ${digest}`);
console.log(`SBOM: ${path.relative(root, sbomPath)}`);
