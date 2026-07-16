import { chmodSync, mkdirSync, writeFileSync } from 'fs';
import { mkdtempSync, rmSync } from 'fs';
import { tmpdir } from 'os';
import { join, resolve } from 'path';
import { spawnSync } from 'child_process';

const scriptPath = resolve(__dirname, '../../../scripts/check-lsp-installed.mjs');

describe('check-lsp-installed script', () => {
  let root: string;
  let codeRoot: string;
  let runtimeRoot: string;
  let registryPath: string;
  let matrixPath: string;
  let fakeBin: string;

  beforeEach(() => {
    root = mkdtempSync(join(tmpdir(), 'openqc-installed-lsp-'));
    codeRoot = join(root, 'code');
    runtimeRoot = join(codeRoot, '.lsp-latest');
    registryPath = join(root, 'registry.ts');
    matrixPath = join(root, 'matrix.json');
    fakeBin = join(root, 'bin');
    mkdirSync(fakeBin, { recursive: true });
    mkdirSync(join(runtimeRoot, '.venv', 'bin'), { recursive: true });

    const registryEntries: string[] = [];
    const readinessEntries: string[] = [];
    const matrixEntries: Record<string, object> = {};
    for (let index = 0; index < 17; index++) {
      const id = `sample-${index}-lsp`;
      const repoName = `sample-${index}`;
      const executable = `sample-${index}-server`;
      const agent = `sample-${index}-tool`;
      registryEntries.push(`{
        id: '${id}',
        repository: 'newtontech/${repoName}',
        executable: '${executable}',
        languageId: 'sample-${index}',
        defaultBranch: 'main',
      }`);
      readinessEntries.push(`'${id}': diagnosticReadiness('${agent}', 'partial')`);
      matrixEntries[id] = { package: id, targetVersion: '1.0.0', channel: 'pypi' };

      const source = join(codeRoot, '.worktrees-lsp-latest', repoName);
      const installed = join(runtimeRoot, repoName);
      mkdirSync(source, { recursive: true });
      mkdirSync(installed, { recursive: true });
      writeFileSync(join(source, 'package.json'), JSON.stringify({ version: '1.0.0' }));
      writeFileSync(join(installed, 'package.json'), JSON.stringify({ version: '1.0.0' }));
      writeFileSync(join(runtimeRoot, '.venv', 'bin', executable), 'server');
      writeFileSync(join(runtimeRoot, '.venv', 'bin', agent), 'agent');
    }
    writeFileSync(
      registryPath,
      `const entries = [${registryEntries.join(',')}];\nconst readiness = {${readinessEntries.join(',')}};\n`
    );
    writeFileSync(
      matrixPath,
      JSON.stringify({ schemaVersion: 'openqc.lsp.release-matrix.v1', entries: matrixEntries })
    );

    const fakeGit = join(fakeBin, 'git');
    writeFileSync(
      fakeGit,
      `#!/bin/sh
path="$2"
command="$3"
if [ "$command" = "rev-parse" ]; then
  if [ -n "$OPENQC_TEST_DRIFT" ] && echo "$path" | grep -q "/.lsp-latest/sample-0$"; then
    printf '%040d\\n' 2
  else
    printf '%040d\\n' 1
  fi
  exit 0
fi
if [ "$command" = "tag" ]; then
  echo "v1.0.0"
  exit 0
fi
exit 1
`
    );
    chmodSync(fakeGit, 0o755);
  });

  afterEach(() => {
    rmSync(root, { recursive: true, force: true });
  });

  function run(extraEnv: NodeJS.ProcessEnv = {}, args = ['--json']) {
    return spawnSync(process.execPath, [scriptPath, ...args], {
      encoding: 'utf8',
      env: {
        ...process.env,
        ...extraEnv,
        PATH: `${fakeBin}:${process.env.PATH}`,
        OPENQC_LSP_CODE_ROOT: codeRoot,
        OPENQC_LSP_RUNTIME_ROOT: runtimeRoot,
        OPENQC_LSP_REGISTRY_PATH: registryPath,
        OPENQC_LSP_RELEASE_MATRIX_PATH: matrixPath,
      },
    });
  }

  it('reports all five provenance fields for a matching installed fleet', () => {
    const result = run();

    expect(result.status).toBe(0);
    const report = JSON.parse(result.stdout);
    expect(report.schemaVersion).toBe('openqc.lsp.runtime-ledger.v1');
    expect(report.summary).toEqual({ total: 17, passing: 17, failing: 0 });
    expect(report.entries[0]).toEqual(
      expect.objectContaining({
        sourceCommit: expect.stringMatching(/^0+1$/),
        installedCommit: expect.stringMatching(/^0+1$/),
        tagVersion: '1.0.0',
        manifestVersion: '1.0.0',
        installedVersion: '1.0.0',
        registryVersion: '1.0.0',
        runtimeParity: 'pass',
      })
    );
  });

  it('fails strict mode when the installed commit drifts from the source', () => {
    const result = run({ OPENQC_TEST_DRIFT: '1' }, ['--json', '--strict']);

    expect(result.status).toBe(1);
    const report = JSON.parse(result.stdout);
    expect(report.summary.failing).toBe(1);
    expect(report.entries[0].runtimeParity).toBe('fail');
    expect(report.entries[0].problems).toContain('installed-commit-drift');
  });
});
