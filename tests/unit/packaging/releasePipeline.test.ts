import { spawnSync } from 'child_process';
import * as fs from 'fs';
import * as os from 'os';
import * as path from 'path';

describe('release pipeline contract', () => {
  const root = path.resolve(__dirname, '../../..');
  const read = (relativePath: string): string =>
    fs.readFileSync(path.join(root, relativePath), 'utf8');
  const packageJson = JSON.parse(read('package.json'));

  it('pins the OpenQC release identity, Node 22, npm, and vsce', () => {
    expect(`${packageJson.publisher}.${packageJson.name}@${packageJson.version}`).toBe(
      'newtontech.openqc@0.0.1'
    );
    expect(packageJson.displayName).toBe('OpenQC - DFT/MD/Quantum Chemistry Suite');
    expect(packageJson.engines.node).toBe('>=22 <23');
    expect(read('.nvmrc').trim()).toBe('22');
    expect(packageJson.packageManager).toBe('npm@10.9.8');
    expect(packageJson.devDependencies['@vscode/vsce']).toBe('3.9.2');
    expect(packageJson.scripts['package:vsix']).not.toContain('npx');
    expect(packageJson.scripts['ci:pr']).toContain('npm run audit:prod');
  });

  it('uses one real Makefile quality gate in CI and release builds', () => {
    const makefile = read('Makefile');
    const ci = read('.github/workflows/ci.yml');
    const release = read('.github/workflows/release.yml');

    expect(makefile).toContain('check: ## Run the canonical pull-request quality gate');
    expect(makefile).toContain('npm run ci:pr');
    expect(makefile).not.toContain('cd core');
    expect(ci).toContain('npx playwright install --with-deps chromium');
    expect(ci).toContain('run: npm run check:release');
    expect(release).toContain('run: npm run check:release');
    expect(ci).toContain('node-version: 22');
    expect(release.match(/node-version: 22/g)).toHaveLength(2);
  });

  it('publishes only after verification and creates the GitHub release last', () => {
    const workflow = read('.github/workflows/release.yml');

    expect(workflow).toContain('build-and-verify:');
    expect(workflow).toContain('publish-marketplace:');
    expect(workflow).toContain('environment: marketplace-production');
    expect(workflow).toContain('needs: build-and-verify');
    expect(workflow).toContain('finalize-github-release:');
    expect(workflow).toContain('needs: publish-marketplace');
    expect(workflow).toContain('vsix/*.vsix.sha256');
    expect(workflow).toContain('vsix/*.sbom.json');
    expect(workflow).toContain('test "${#vsix_files[@]}" -eq 1');
  });

  it('guards clean origin/master tags and rejects VSIX development payloads', () => {
    const guard = read('scripts/release-guard.mjs');
    const verifier = read('scripts/verify-vsix.mjs');
    const vscodeignore = read('.vscodeignore');

    expect(guard).toContain("['status', '--porcelain', '--untracked-files=normal']");
    expect(guard).toContain("['rev-parse', 'refs/remotes/origin/master']");
    expect(guard).toContain('tagCommit !== head');
    expect(verifier).toContain('VSIX contains forbidden files');
    expect(verifier).toContain('extension/media/vendor/3dmol/3Dmol-min.js');
    expect(verifier).toContain('extension/media/vendor/ngl/ngl.js');
    expect(verifier).toContain('extension/media/vendor/plotly.js-dist-min/plotly.min.js');
    expect(verifier).not.toContain('extension/node_modules/3dmol');
    expect(vscodeignore).toContain('core/**');
    expect(vscodeignore).toContain('tests/**');
    expect(vscodeignore).not.toMatch(/^!node_modules\/(?:3dmol|ngl|plotly\.js-dist-min)\/\*\*$/m);
  });

  it('retires the duplicate Python core after preserving the maintained workflow implementation', () => {
    expect(fs.existsSync(path.join(root, 'core/pyproject.toml'))).toBe(false);
    expect(fs.existsSync(path.join(root, 'core/openqc/workflow.py'))).toBe(false);
    expect(fs.existsSync(path.join(root, 'src/utils/migration/mdWorkflow.ts'))).toBe(true);
    expect(fs.existsSync(path.join(root, 'tests/unit/migration/mdWorkflow.test.ts'))).toBe(true);
  });

  it('enforces clean, synchronized, version-matched release tags', () => {
    const fixtureRoot = fs.mkdtempSync(path.join(os.tmpdir(), 'openqc-release-guard-'));
    const repository = path.join(fixtureRoot, 'repository');
    const remote = path.join(fixtureRoot, 'remote.git');
    const gitEnv = { ...process.env };
    for (const variable of ['GIT_DIR', 'GIT_WORK_TREE', 'GIT_INDEX_FILE', 'GIT_PREFIX']) {
      delete gitEnv[variable];
    }
    fs.mkdirSync(repository);

    const git = (cwd: string, args: string[]) => {
      const result = spawnSync('git', args, { cwd, encoding: 'utf8', env: gitEnv });
      if (result.status !== 0) {
        throw new Error(result.stderr || result.stdout);
      }
    };
    const guard = (...args: string[]) =>
      spawnSync(process.execPath, [path.join(root, 'scripts/release-guard.mjs'), ...args], {
        cwd: repository,
        encoding: 'utf8',
        env: { ...gitEnv, OPENQC_RELEASE_ROOT: repository },
      });

    try {
      fs.writeFileSync(path.join(repository, 'package.json'), '{"version":"0.0.1"}\n');
      git(repository, ['init', '-b', 'master']);
      git(repository, ['config', 'user.name', 'OpenQC Test']);
      git(repository, ['config', 'user.email', 'openqc-test@example.invalid']);
      git(repository, ['add', 'package.json']);
      git(repository, ['commit', '-m', 'fixture']);
      git(fixtureRoot, ['init', '--bare', remote]);
      git(repository, ['remote', 'add', 'origin', remote]);
      git(repository, ['push', '-u', 'origin', 'master']);

      expect(guard('--prepare', '--no-fetch').status).toBe(0);

      fs.writeFileSync(path.join(repository, 'dirty.txt'), 'dirty\n');
      const dirty = guard('--prepare', '--no-fetch');
      expect(dirty.status).toBe(1);
      expect(dirty.stderr).toContain('working tree is not clean');
      fs.rmSync(path.join(repository, 'dirty.txt'));

      git(repository, ['tag', '-a', 'v0.0.1', '-m', 'Release v0.0.1']);
      expect(guard('--tag', 'v0.0.1', '--no-fetch').status).toBe(0);

      fs.writeFileSync(path.join(repository, 'ahead.txt'), 'ahead\n');
      git(repository, ['add', 'ahead.txt']);
      git(repository, ['commit', '-m', 'ahead']);
      const ahead = guard('--tag', 'v0.0.1', '--no-fetch');
      expect(ahead.status).toBe(1);
      expect(ahead.stderr).toContain('is not synchronized with origin/master');

      expect(guard('--tag', 'v9.9.9', '--no-fetch').stderr).toContain(
        'must match package.json version as v0.0.1'
      );
    } finally {
      fs.rmSync(fixtureRoot, { recursive: true, force: true });
    }
  });
});
