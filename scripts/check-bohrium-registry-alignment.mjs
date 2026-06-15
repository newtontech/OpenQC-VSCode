#!/usr/bin/env node

/**
 * OpenQC ↔ Bohrium LSP registry alignment check.
 *
 * Compares bundled OpenQC registry ids/agent CLIs with the Bohrium skill
 * backend registry. Surfaces missing or excess entries without conflating
 * OpenQC editor integration readiness with Bohrium skill-router availability.
 *
 * Usage:
 *   node scripts/check-bohrium-registry-alignment.mjs
 *   node scripts/check-bohrium-registry-alignment.mjs --json
 *   node scripts/check-bohrium-registry-alignment.mjs --strict
 *   node scripts/check-bohrium-registry-alignment.mjs --bohrium-registry /path/to/lsp_backends.yaml
 *
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/194
 */

import { existsSync, readFileSync } from 'fs';
import { dirname, join, resolve } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '..');
const codeRoot = resolve(repoRoot, '..');

const jsonMode = process.argv.includes('--json');
const strictMode = process.argv.includes('--strict');
const registryArgIndex = process.argv.indexOf('--bohrium-registry');
const defaultBohriumRegistry = join(
  codeRoot,
  'bohrium_skills/bohrium_skills/lsp-skills/references/lsp_backends.yaml'
);
const bohriumRegistryPath =
  registryArgIndex >= 0 ? resolve(process.argv[registryArgIndex + 1]) : defaultBohriumRegistry;

const openqcRegistryPath = join(repoRoot, 'src/lsp/registry.ts');
const openqcRegistry = readFileSync(openqcRegistryPath, 'utf8');

const openqcEntries = [
  ...openqcRegistry.matchAll(/id: '([^']+)'[\s\S]*?languageId: '([^']+)'/g),
].map(match => ({
  id: match[1],
  languageId: match[2],
}));

const openqcReadiness = new Map(
  [...openqcRegistry.matchAll(/'([^']+)': diagnosticReadiness\('([^']+)', '([^']+)'/g)].map(
    match => [
      match[1],
      {
        agentCli: match[2],
        closedLoop: match[3],
      },
    ]
  )
);

if (!existsSync(bohriumRegistryPath)) {
  const message = `Bohrium registry not found: ${bohriumRegistryPath}`;
  if (jsonMode) {
    console.log(JSON.stringify({ ok: false, error: message }, null, 2));
  } else {
    console.error(message);
  }
  process.exit(strictMode ? 1 : 0);
}

const bohriumRegistry = JSON.parse(readFileSync(bohriumRegistryPath, 'utf8'));
const bohriumBackends = (bohriumRegistry.backends ?? []).map(backend => ({
  id: backend.id,
  software: backend.software,
  agentCli: backend.agent_cli ?? backend.executable,
  maturity: backend.maturity,
  blocking: Boolean(backend.blocking),
}));

const openqcIds = new Set(openqcEntries.map(entry => entry.id));
const bohriumIds = new Set(bohriumBackends.map(backend => backend.id));

const missingInBohrium = [...openqcIds].filter(id => !bohriumIds.has(id)).sort();
const excessInBohrium = [...bohriumIds].filter(id => !openqcIds.has(id)).sort();

const agentCliMismatches = [];
for (const entry of openqcEntries) {
  const backend = bohriumBackends.find(item => item.id === entry.id);
  if (!backend) {
    continue;
  }
  const expectedAgentCli = openqcReadiness.get(entry.id)?.agentCli;
  if (expectedAgentCli && backend.agentCli !== expectedAgentCli) {
    agentCliMismatches.push({
      id: entry.id,
      openqcAgentCli: expectedAgentCli,
      bohriumAgentCli: backend.agentCli,
    });
  }
}

const aligned =
  missingInBohrium.length === 0 && excessInBohrium.length === 0 && agentCliMismatches.length === 0;

const report = {
  ok: aligned,
  generatedAt: new Date().toISOString(),
  surfaces: {
    openqc: {
      role: 'editor-integration',
      source: openqcRegistryPath,
      backendCount: openqcEntries.length,
      description: 'Bundled VS Code language contributions and local LSP startup registry.',
    },
    bohrium: {
      role: 'skill-router',
      source: bohriumRegistryPath,
      backendCount: bohriumBackends.length,
      description:
        'Bohrium LSP skill backend registry used by lsp_agent_router.py and lsp_skill_probe.py.',
    },
  },
  summary: {
    openqcCount: openqcEntries.length,
    bohriumCount: bohriumBackends.length,
    missingInBohrium,
    excessInBohrium,
    agentCliMismatches,
  },
  note: 'OpenQC editor integration and Bohrium skill routing share backend ids and agent CLI names, but maturity/blocking in Bohrium may differ from OpenQC stability until fleet closeout.',
  backends: openqcEntries.map(entry => {
    const readiness = openqcReadiness.get(entry.id);
    const backend = bohriumBackends.find(item => item.id === entry.id);
    return {
      id: entry.id,
      languageId: entry.languageId,
      openqcAgentCli: readiness?.agentCli ?? null,
      openqcClosedLoop: readiness?.closedLoop ?? null,
      bohriumPresent: Boolean(backend),
      bohriumAgentCli: backend?.agentCli ?? null,
      bohriumMaturity: backend?.maturity ?? null,
      bohriumBlocking: backend?.blocking ?? null,
    };
  }),
};

if (jsonMode) {
  console.log(JSON.stringify(report, null, 2));
} else {
  console.log('# OpenQC ↔ Bohrium LSP Registry Alignment');
  console.log('');
  console.log(
    `OpenQC editor registry: ${report.surfaces.openqc.backendCount} backends (${openqcRegistryPath})`
  );
  console.log(
    `Bohrium skill registry: ${report.surfaces.bohrium.backendCount} backends (${bohriumRegistryPath})`
  );
  console.log('');
  if (aligned) {
    console.log('Status: aligned');
  } else {
    console.log('Status: drift detected');
    if (missingInBohrium.length > 0) {
      console.log(`- Missing in Bohrium: ${missingInBohrium.join(', ')}`);
    }
    if (excessInBohrium.length > 0) {
      console.log(`- Excess in Bohrium: ${excessInBohrium.join(', ')}`);
    }
    if (agentCliMismatches.length > 0) {
      console.log(`- Agent CLI mismatches: ${agentCliMismatches.map(item => item.id).join(', ')}`);
    }
  }
  console.log('');
  console.log(report.note);
}

if (strictMode && !aligned) {
  process.exit(1);
}
