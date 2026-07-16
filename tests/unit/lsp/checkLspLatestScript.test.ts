import { spawnSync } from 'child_process';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

const REPO_ROOT = path.resolve(__dirname, '../../..');
const SCRIPT_PATH = path.join(REPO_ROOT, 'scripts', 'check-lsp-latest.mjs');

describe('check-lsp-latest script remote fallback', () => {
  let tempDir: string;

  beforeEach(() => {
    tempDir = fs.mkdtempSync(path.join(os.tmpdir(), 'openqc-lsp-latest-'));
  });

  afterEach(() => {
    fs.rmSync(tempDir, { recursive: true, force: true });
  });

  it('tries SSH remote HEAD lookup when HTTPS ls-remote times out', () => {
    const head = 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';
    const harness = createHarness(tempDir, {
      localHead: head,
      httpsError: 'fatal: unable to access: Operation timed out',
      sshHead: head,
    });

    const result = runLatestCheck(harness);

    expect(result.status).toBe(0);
    expect(result.stdout).toContain('fake-lsp');
    expect(result.stdout).toContain('latest');
    expect(result.stdout).not.toContain('remote-unavailable');
    expect(readGitCalls(harness)).toContain('ls-remote https://github.com/newtontech/fake-lsp.git');
    expect(readGitCalls(harness)).toContain('ls-remote git@github.com:newtontech/fake-lsp.git');
  });

  it('does not try SSH when HTTPS succeeds but local checkout is stale', () => {
    const localHead = 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa';
    const remoteHead = 'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb';
    const harness = createHarness(tempDir, { localHead, httpsHead: remoteHead });

    const result = runLatestCheck(harness);

    expect(result.status).toBe(1);
    expect(result.stdout).toContain('fake-lsp');
    expect(result.stdout).toContain('not-at-remote-head');
    expect(result.stdout).not.toContain('remote-unavailable');
    expect(readGitCalls(harness)).toContain('ls-remote https://github.com/newtontech/fake-lsp.git');
    expect(readGitCalls(harness)).not.toContain('ls-remote git@github.com:newtontech/fake-lsp.git');
  });

  it('does not try SSH for non-network HTTPS ls-remote errors', () => {
    const harness = createHarness(tempDir, {
      localHead: 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa',
      httpsError: 'fatal: repository not found',
    });

    const result = runLatestCheck(harness);

    expect(result.status).toBe(1);
    expect(result.stdout).toContain('remote-unavailable');
    expect(result.stdout).toContain('fake-lsp');
    expect(readGitCalls(harness)).toContain('ls-remote https://github.com/newtontech/fake-lsp.git');
    expect(readGitCalls(harness)).not.toContain('ls-remote git@github.com:newtontech/fake-lsp.git');
  });
});

interface HarnessOptions {
  localHead: string;
  httpsHead?: string;
  httpsError?: string;
  sshHead?: string;
}

interface Harness {
  registryPath: string;
  codeRoot: string;
  binDir: string;
  gitLog: string;
}

function createHarness(tempDir: string, options: HarnessOptions): Harness {
  const registryPath = path.join(tempDir, 'registry.ts');
  const codeRoot = path.join(tempDir, 'code');
  const repoRoot = path.join(codeRoot, 'fake-lsp');
  const binDir = path.join(tempDir, 'bin');
  const gitLog = path.join(tempDir, 'git.log');

  fs.mkdirSync(path.join(repoRoot, '.git'), { recursive: true });
  fs.mkdirSync(binDir, { recursive: true });
  fs.writeFileSync(
    registryPath,
    [
      'export const LSP_SERVERS = [{',
      "  id: 'fake-lsp',",
      "  repository: 'newtontech/fake-lsp',",
      "  languageId: 'fake',",
      "  defaultBranch: 'main',",
      '}];',
    ].join('\n')
  );
  writeFakeGit(binDir, gitLog, options);

  return { registryPath, codeRoot, binDir, gitLog };
}

function runLatestCheck(harness: Harness): ReturnType<typeof spawnSync> {
  return spawnSync(process.execPath, [SCRIPT_PATH, '--fail-on-drift'], {
    cwd: REPO_ROOT,
    encoding: 'utf8',
    env: {
      ...process.env,
      PATH: `${harness.binDir}${path.delimiter}${process.env.PATH ?? ''}`,
      OPENQC_LSP_REGISTRY_PATH: harness.registryPath,
      OPENQC_LSP_CODE_ROOT: harness.codeRoot,
    },
  });
}

function readGitCalls(harness: Harness): string {
  return fs.readFileSync(harness.gitLog, 'utf8');
}

function shellSingleQuote(value: string): string {
  return `'${value.replace(/'/g, "'\\''")}'`;
}

function writeFakeGit(binDir: string, logPath: string, options: HarnessOptions): void {
  const fakeGit = path.join(binDir, 'git');
  const httpsBranch = options.httpsHead
    ? `    https://*) printf '${options.httpsHead}\\trefs/heads/main\\n'; exit 0 ;;`
    : `    https://*) echo ${shellSingleQuote(options.httpsError ?? 'fatal: unable to access')} >&2; exit 128 ;;`;
  const sshBranch = options.sshHead
    ? `    git@github.com:*) printf '${options.sshHead}\\trefs/heads/main\\n'; exit 0 ;;`
    : '    git@github.com:*) echo "unexpected SSH fallback" >&2; exit 2 ;;';

  fs.writeFileSync(
    fakeGit,
    [
      '#!/bin/sh',
      `printf '%s\\n' "$*" >> "${logPath}"`,
      'if [ "$1" = "ls-remote" ]; then',
      '  case "$2" in',
      httpsBranch,
      sshBranch,
      '  esac',
      'fi',
      'if [ "$1" = "-C" ]; then',
      '  if [ "$3" = "rev-parse" ]; then',
      `    echo "${options.localHead}"; exit 0`,
      '  fi',
      '  if [ "$3" = "status" ]; then',
      '    exit 0',
      '  fi',
      'fi',
      'echo "unexpected git invocation: $*" >&2',
      'exit 2',
    ].join('\n'),
    { mode: 0o755 }
  );
}
