#!/usr/bin/env node

import { execFileSync, spawnSync } from 'child_process';
import {
  chmodSync,
  existsSync,
  mkdirSync,
  readFileSync,
  rmSync,
  symlinkSync,
  writeFileSync,
} from 'fs';
import { dirname, join, resolve } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '..');
const codeRoot = resolve(repoRoot, '..');
const registryPath = join(repoRoot, 'src/lsp/registry.ts');
const managedRoot = join(codeRoot, '.lsp-latest');
const latestWorktreeRoot = join(codeRoot, '.worktrees-lsp-latest');
const buildRoot = join(managedRoot, '.build');
const venvDir = join(managedRoot, '.venv');
const localBin = join(process.env.HOME || '', '.local', 'bin');
const cargoBin = join(process.env.HOME || '', '.cargo', 'bin');
const dryRun = process.argv.includes('--dry-run');
const skipInstall = process.argv.includes('--skip-install');

const registry = readFileSync(registryPath, 'utf8');
const entries = extractRegistryEntryBlocks(registry).map(block => ({
  id: readStringField(block, 'id'),
  repository: readStringField(block, 'repository'),
  executable: readStringField(block, 'executable'),
  languageId: readStringField(block, 'languageId'),
  defaultBranch: readStringField(block, 'defaultBranch'),
  repoName: readStringField(block, 'repoName'),
  localLaunchKind: readNestedStringField(block, 'localLaunch', 'kind'),
  nodeScriptPath: readNestedStringField(block, 'localLaunch', 'scriptPath'),
  cargoBinaryName: readNestedStringField(block, 'localLaunch', 'binaryName'),
  agentCli: readAgentCli(registry, readStringField(block, 'id')),
}));

if (entries.length !== 17) {
  throw new Error(`Expected 17 LSP registry entries, found ${entries.length}`);
}

mkdirSync(managedRoot, { recursive: true });
mkdirSync(latestWorktreeRoot, { recursive: true });
mkdirSync(buildRoot, { recursive: true });
mkdirSync(localBin, { recursive: true });
mkdirSync(cargoBin, { recursive: true });

const rows = [];

if (!skipInstall) {
  ensurePythonVenv();
}

for (const entry of entries) {
  const managedPath = join(managedRoot, entry.repoName);
  const latestPath = join(latestWorktreeRoot, entry.repoName);
  const remoteUrl = `https://github.com/${entry.repository}.git`;

  refreshCheckout(remoteUrl, entry.defaultBranch, managedPath);
  refreshCheckout(remoteUrl, entry.defaultBranch, latestPath);

  const head = git(['-C', managedPath, 'rev-parse', '--short=12', 'HEAD']);
  const installed = skipInstall ? [] : installEntry(entry, managedPath);
  const commandPath = resolveCommand(entry.executable) || '-';
  const agentPath = entry.agentCli ? resolveCommand(entry.agentCli) || '-' : '-';
  const probe = commandPath !== '-' || agentPath !== '-' ? 'ok' : 'missing';

  rows.push({
    id: entry.id,
    language: entry.languageId,
    branch: entry.defaultBranch,
    head,
    kind: entry.localLaunchKind || 'command',
    command: commandPath,
    agentCli: agentPath,
    installed: installed.join(', ') || '-',
    probe,
  });
}

console.table(rows);

const missing = rows.filter(row => row.probe !== 'ok');
if (missing.length > 0) {
  console.error(`Missing installed commands for ${missing.length} LSP backend(s).`);
  process.exitCode = 1;
}

function installEntry(entry, checkoutPath) {
  if (entry.localLaunchKind === 'nodeScript') {
    return installNodeEntry(entry, checkoutPath);
  }

  if (entry.localLaunchKind === 'cargoBinary') {
    return installCargoEntry(entry, checkoutPath);
  }

  return installPythonEntry(entry, checkoutPath);
}

function installPythonEntry(entry, checkoutPath) {
  run(venvPython(), ['-m', 'pip', 'install', '--upgrade', 'pip', 'wheel', 'setuptools']);
  const installPath = preparePythonInstallPath(entry, checkoutPath);
  run(venvPython(), ['-m', 'pip', 'install', '--upgrade', installPath]);

  const linked = [];
  for (const command of unique([entry.executable, entry.agentCli])) {
    const source = join(venvDir, 'bin', command);
    if (existsSync(source)) {
      linkCommand(source, join(localBin, command));
      linked.push(command);
    }
  }
  return linked;
}

function preparePythonInstallPath(entry, checkoutPath) {
  if (entry.id !== 'cp2k-lsp-enhanced') {
    return checkoutPath;
  }

  const buildPath = join(buildRoot, entry.repoName);
  rmSync(buildPath, { recursive: true, force: true });
  run('rsync', ['-a', '--exclude', '.git', `${checkoutPath}/`, `${buildPath}/`]);
  removeDuplicatePoetryScripts(join(buildPath, 'pyproject.toml'));
  return buildPath;
}

function removeDuplicatePoetryScripts(pyprojectPath) {
  if (!existsSync(pyprojectPath)) {
    return;
  }

  const lines = readFileSync(pyprojectPath, 'utf8').split('\n');
  const seen = new Set();
  let inScripts = false;
  const cleaned = [];

  for (const line of lines) {
    const section = line.match(/^\[([^\]]+)\]\s*$/);
    if (section) {
      inScripts = section[1] === 'tool.poetry.scripts';
      cleaned.push(line);
      continue;
    }

    if (inScripts) {
      const key = line.match(/^([A-Za-z0-9_.-]+)\s*=/)?.[1];
      if (key) {
        if (seen.has(key)) {
          continue;
        }
        seen.add(key);
      }
    }

    cleaned.push(line);
  }

  writeFileSync(pyprojectPath, cleaned.join('\n'));
}

function installNodeEntry(entry, checkoutPath) {
  if (existsSync(join(checkoutPath, 'package-lock.json'))) {
    run('npm', ['ci'], { cwd: checkoutPath });
  } else {
    run('npm', ['install'], { cwd: checkoutPath });
  }

  const pkgPath = join(checkoutPath, 'package.json');
  const pkg = existsSync(pkgPath) ? JSON.parse(readFileSync(pkgPath, 'utf8')) : {};
  if (pkg.scripts?.compile) {
    run('npm', ['run', 'compile'], { cwd: checkoutPath });
  } else if (pkg.scripts?.build) {
    run('npm', ['run', 'build'], { cwd: checkoutPath });
  }

  const linked = [];
  const binEntries = normalizePackageBins(pkg.bin);
  for (const [name, relativePath] of Object.entries(binEntries)) {
    const source = join(checkoutPath, relativePath);
    if (existsSync(source)) {
      linkCommand(source, join(localBin, name));
      linked.push(name);
    }
  }

  if (entry.executable && !linked.includes(entry.executable)) {
    const target = resolveNodeScript(entry, checkoutPath);
    if (target) {
      writeNodeWrapper(join(localBin, entry.executable), target);
      linked.push(entry.executable);
    }
  }

  if (entry.agentCli && !linked.includes(entry.agentCli)) {
    const target =
      resolveNodeScript({ ...entry, nodeScriptPath: 'dist/tool.js' }, checkoutPath) ||
      resolveNodeScript({ ...entry, nodeScriptPath: 'server/out/tool.js' }, checkoutPath) ||
      resolveNodeScript(entry, checkoutPath);
    if (target) {
      writeNodeWrapper(join(localBin, entry.agentCli), target);
      linked.push(entry.agentCli);
    }
  }

  return linked;
}

function installCargoEntry(entry, checkoutPath) {
  run('cargo', ['install', '--locked', '--path', checkoutPath, '--force']);
  return unique([entry.executable, entry.agentCli]).filter(command => resolveCommand(command));
}

function refreshCheckout(remoteUrl, branch, checkoutPath) {
  if (!existsSync(join(checkoutPath, '.git'))) {
    if (existsSync(checkoutPath)) {
      rmSync(checkoutPath, { recursive: true, force: true });
    }
    run('git', ['clone', '--branch', branch, '--single-branch', remoteUrl, checkoutPath]);
    return;
  }

  run('git', ['-C', checkoutPath, 'fetch', 'origin', branch, '--prune']);
  run('git', ['-C', checkoutPath, 'reset', '--hard', 'HEAD']);
  run('git', ['-C', checkoutPath, 'clean', '-fdx']);
  run('git', ['-C', checkoutPath, 'checkout', '-B', branch, `origin/${branch}`]);
  run('git', ['-C', checkoutPath, 'reset', '--hard', `origin/${branch}`]);
  run('git', ['-C', checkoutPath, 'clean', '-fdx']);
}

function ensurePythonVenv() {
  if (!existsSync(join(venvDir, 'bin', 'python'))) {
    run('python3', ['-m', 'venv', venvDir]);
  }
}

function venvPython() {
  return join(venvDir, 'bin', 'python');
}

function run(command, args, options = {}) {
  const display = [command, ...args].join(' ');
  console.log(`$ ${display}`);
  if (dryRun) {
    return '';
  }

  const result = spawnSync(command, args, {
    cwd: options.cwd,
    env: {
      ...process.env,
      PATH: `${localBin}:${cargoBin}:${process.env.PATH || ''}`,
    },
    encoding: 'utf8',
    stdio: options.capture ? ['ignore', 'pipe', 'pipe'] : 'inherit',
  });

  if (result.status !== 0) {
    const stderr = result.stderr ? `\n${result.stderr}` : '';
    throw new Error(`Command failed (${result.status}): ${display}${stderr}`);
  }

  return result.stdout?.trim() || '';
}

function git(args) {
  return execFileSync('git', args, { encoding: 'utf8' }).trim();
}

function linkCommand(source, target) {
  chmodSync(source, 0o755);
  rmSync(target, { force: true });
  symlinkSync(source, target);
}

function writeNodeWrapper(target, scriptPath) {
  rmSync(target, { force: true });
  writeFileSync(target, `#!/usr/bin/env node\nimport '${scriptPath.replaceAll("'", '%27')}';\n`);
  run('chmod', ['755', target]);
}

function resolveNodeScript(entry, checkoutPath) {
  const candidates = unique([
    entry.nodeScriptPath,
    'server/out/server.js',
    'out/server.js',
    'dist/server.js',
    'dist/index.js',
    'index.js',
  ]).filter(Boolean);

  for (const candidate of candidates) {
    const absolute = join(checkoutPath, candidate);
    if (existsSync(absolute)) {
      return absolute;
    }
  }
  return undefined;
}

function resolveCommand(command) {
  if (!command) {
    return '';
  }
  try {
    return execFileSync('which', [command], {
      encoding: 'utf8',
      env: {
        ...process.env,
        PATH: `${localBin}:${cargoBin}:${process.env.PATH || ''}`,
      },
    }).trim();
  } catch {
    return '';
  }
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

function readStringField(block, field) {
  const match = block.match(new RegExp(`${field}: '([^']+)'`));
  return match?.[1] ?? '';
}

function readNestedStringField(block, parent, field) {
  const match = block.match(new RegExp(`${parent}: \\{([\\s\\S]*?)\\n\\s*\\}`));
  return match ? readStringField(match[1], field) : '';
}

function readAgentCli(text, id) {
  const match = text.match(new RegExp(`'${escapeRegExp(id)}': \\{[\\s\\S]*?agentCli: '([^']+)'`));
  return match?.[1] ?? '';
}

function normalizePackageBins(bin) {
  if (!bin) {
    return {};
  }
  if (typeof bin === 'string') {
    return { cli: bin };
  }
  return bin;
}

function unique(values) {
  return [...new Set(values.filter(Boolean))];
}

function escapeRegExp(value) {
  return value.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}
