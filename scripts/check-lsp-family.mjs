#!/usr/bin/env node

import { execFileSync } from 'child_process';
import { existsSync, readFileSync } from 'fs';
import { dirname, join, relative, resolve } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '..');
const codeRoot = resolve(repoRoot, '..');
const registryPath = join(repoRoot, 'src/lsp/registry.ts');
const packagePath = join(repoRoot, 'package.json');
const asJson = process.argv.includes('--json');
const failOnWarn = process.argv.includes('--fail-on-warn');

const requiredOperations = [
  'capabilities',
  'check',
  'context',
  'complete',
  'hover',
  'symbols',
  'fix',
];
const requiredCapabilities = [
  'diagnostics',
  'rich-diagnostics',
  'completion',
  'hover',
  'symbols',
  'fix-preview',
  'llm-wiki',
  'openqc-context',
];

const registry = readFileSync(registryPath, 'utf8');
const pkg = JSON.parse(readFileSync(packagePath, 'utf8'));
const entries = extractRegistryEntryBlocks(registry).map(block => ({
  id: readStringField(block, 'id'),
  repository: readStringField(block, 'repository'),
  executable: readStringField(block, 'executable'),
  languageId: readStringField(block, 'languageId'),
  fileExtensions: parseArrayField(block, 'fileExtensions'),
  fileNames: parseArrayField(block, 'fileNames'),
  defaultBranch: readStringField(block, 'defaultBranch'),
  repoName: readStringField(block, 'repoName'),
}));

if (entries.length !== 17) {
  throw new Error(`Expected 17 registry entries, found ${entries.length}`);
}

const languages = pkg.contributes?.languages ?? [];
const properties = pkg.contributes?.configuration?.properties ?? {};

const rows = entries.map(entry => checkEntry(entry));
const failed = rows.filter(row => row.status === 'fail');
const warned = rows.filter(row => row.status === 'warn');

if (asJson) {
  console.log(
    JSON.stringify(
      {
        generatedAt: new Date().toISOString(),
        repoRoot,
        codeRoot,
        total: rows.length,
        failed: failed.length,
        warned: warned.length,
        rows,
      },
      null,
      2
    )
  );
} else {
  console.table(
    rows.map(row => ({
      id: row.id,
      language: row.languageId,
      repo: row.repoExists ? 'yes' : 'no',
      manifest: row.manifestExists ? 'yes' : 'no',
      cli: row.agentCli || '-',
      commit: row.commit || '-',
      valid: row.smokeFixtures.valid ? 'yes' : 'no',
      invalid: row.smokeFixtures.invalid ? 'yes' : 'no',
      logs: row.smokeFixtures.logs ? 'yes' : 'no',
      status: row.status,
    }))
  );
  for (const row of rows) {
    if (row.issues.length > 0) {
      console.log(`\\n${row.id}`);
      for (const issue of row.issues) {
        console.log(`  - ${issue}`);
      }
    }
  }
}

if (failed.length > 0 || (failOnWarn && warned.length > 0)) {
  process.exitCode = 1;
}

function checkEntry(entry) {
  const repoPath = join(codeRoot, entry.repoName);
  const manifestPath = join(repoPath, 'lsp-capabilities.json');
  const issues = [];
  const language = languages.find(item => item.id === entry.languageId);
  if (!language) {
    issues.push(`package.json missing language ${entry.languageId}`);
  }
  for (const key of [
    `openqc.lsp.${entry.languageId}.enabled`,
    `openqc.lsp.${entry.languageId}.path`,
    `openqc.lsp.${entry.languageId}.command`,
    `openqc.lsp.${entry.languageId}.args`,
    `openqc.lsp.${entry.languageId}.env`,
  ]) {
    if (!(key in properties)) {
      issues.push(`package.json missing setting ${key}`);
    }
  }

  let manifest = null;
  if (!existsSync(manifestPath)) {
    issues.push(`missing ${relative(repoRoot, manifestPath)}`);
  } else {
    try {
      manifest = JSON.parse(readFileSync(manifestPath, 'utf8'));
      validateManifest(entry, manifest, repoPath, issues);
    } catch (error) {
      issues.push(`invalid manifest JSON: ${error.message}`);
    }
  }

  const commit = git(['-C', repoPath, 'rev-parse', '--short=12', 'HEAD']);
  const dirty = git(['-C', repoPath, 'status', '--short']);
  if (!existsSync(repoPath)) {
    issues.push(`missing sibling repo ${repoPath}`);
  }

  const hardFailures = issues.filter(
    issue =>
      issue.startsWith('missing sibling repo') ||
      issue.startsWith('missing ') ||
      issue.startsWith('invalid manifest JSON') ||
      issue.includes('must ') ||
      issue.includes('!=') ||
      issue.includes('missing operation')
  );

  return {
    id: entry.id,
    languageId: entry.languageId,
    repoPath,
    repoExists: existsSync(repoPath),
    commit,
    dirty: dirty ? 'yes' : 'no',
    manifestPath,
    manifestExists: existsSync(manifestPath),
    agentCli: manifest?.agentCli?.command,
    capabilities: manifest?.capabilities ?? [],
    smokeFixtures: {
      valid: hasAny(repoPath, manifest?.fixturePaths?.valid ?? []),
      invalid: hasAny(repoPath, manifest?.fixturePaths?.invalid ?? []),
      logs: hasAny(repoPath, manifest?.fixturePaths?.logs ?? []),
    },
    issues,
    status: hardFailures.length > 0 ? 'fail' : issues.length > 0 ? 'warn' : 'pass',
  };
}

function validateManifest(entry, manifest, repoPath, issues) {
  if (manifest.schema !== 'OpenQCLspCapabilities')
    issues.push('schema must be OpenQCLspCapabilities');
  if (manifest.version !== 1) issues.push('version must be 1');
  if (manifest.id !== entry.id) issues.push(`manifest id ${manifest.id} != ${entry.id}`);
  if (manifest.languageId !== entry.languageId)
    issues.push(`manifest languageId ${manifest.languageId} != ${entry.languageId}`);
  if (manifest.executable !== entry.executable)
    issues.push(`manifest executable ${manifest.executable} != ${entry.executable}`);
  if (manifest.defaultBranch !== entry.defaultBranch)
    issues.push(`manifest defaultBranch ${manifest.defaultBranch} != ${entry.defaultBranch}`);
  if (!Array.isArray(manifest.filePatterns) || manifest.filePatterns.length === 0)
    issues.push('filePatterns must be non-empty');
  if (!['blocking', 'warning-only'].includes(manifest.blockingPolicy?.mode))
    issues.push('blockingPolicy.mode must be blocking or warning-only');
  if (!Array.isArray(manifest.capabilities) || manifest.capabilities.length === 0)
    issues.push('capabilities must be non-empty');
  for (const capability of requiredCapabilities) {
    if (!manifest.capabilities?.includes(capability))
      issues.push(`missing capability ${capability}`);
  }
  if (!manifest.agentCli?.command) issues.push('agentCli.command must be set');
  for (const operation of requiredOperations) {
    if (!manifest.agentCli?.operations?.includes(operation))
      issues.push(`missing operation ${operation}`);
  }
  if (manifest.agentCli?.jsonFormat !== true) issues.push('agentCli.jsonFormat must be true');
  if (manifest.agentCli?.failOnBlocking !== true)
    issues.push('agentCli.failOnBlocking must be true');
  for (const key of ['diagnosticSchema']) {
    if (!manifest[key] || !existsSync(join(repoPath, manifest[key])))
      issues.push(`missing manifest path ${manifest[key] || key}`);
  }
  for (const key of ['plan', 'rawAssets', 'entities', 'concepts', 'synthesis', 'index', 'log']) {
    const value = manifest.wikiPaths?.[key];
    if (!value || !existsSync(join(repoPath, value)))
      issues.push(`missing wiki path ${value || key}`);
  }
  if (manifest.openqc?.registryId !== entry.id)
    issues.push('openqc.registryId must match registry id');
  if (manifest.openqc?.contextContract !== 'DSLAuthoringContext')
    issues.push('openqc.contextContract must be DSLAuthoringContext');
  if (manifest.openqc?.diagnosticEnvelope !== 'DiagnosticEnvelope/v1')
    issues.push('openqc.diagnosticEnvelope must be DiagnosticEnvelope/v1');
}

function parseStringArray(text) {
  return [...text.matchAll(/'([^']+)'/g)].map(match => match[1]);
}

function extractRegistryEntryBlocks(text) {
  const start = text.indexOf('const BUNDLED_LSP_SERVERS');
  const arrayStart = text.indexOf('[', start);
  const arrayEnd = text.indexOf('] as const', arrayStart);
  if (start === -1 || arrayStart === -1 || arrayEnd === -1) {
    return [];
  }
  const body = text.slice(arrayStart + 1, arrayEnd);
  const blocks = [];
  let depth = 0;
  let blockStart = -1;
  for (let index = 0; index < body.length; index += 1) {
    const char = body[index];
    if (char === '{') {
      if (depth === 0) {
        blockStart = index;
      }
      depth += 1;
    } else if (char === '}') {
      depth -= 1;
      if (depth === 0 && blockStart !== -1) {
        blocks.push(body.slice(blockStart, index + 1));
        blockStart = -1;
      }
    }
  }
  return blocks.filter(block => block.includes("id: '"));
}

function parseArrayField(block, field) {
  const match = block.match(new RegExp(`${field}: \\[([\\s\\S]*?)\\]`));
  return match ? parseStringArray(match[1]) : [];
}

function readStringField(block, field) {
  const match = block.match(new RegExp(`${field}: '([^']+)'`));
  return match?.[1] ?? '';
}

function hasAny(repoPath, candidates) {
  return candidates.some(candidate => existsSync(join(repoPath, candidate)));
}

function git(args) {
  try {
    return execFileSync('git', args, {
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'ignore'],
    }).trim();
  } catch {
    return '';
  }
}
