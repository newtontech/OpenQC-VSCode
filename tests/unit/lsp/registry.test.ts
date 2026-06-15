import * as fs from 'fs';
import * as path from 'path';
import packageJson from '../../../package.json';
import {
  getBundledLspServerCount,
  getLspDiagnosticReadiness,
  getLspServerByLanguageId,
  getLspServerBySoftwareName,
  listBundledLspServers,
} from '../../../src/lsp/registry';
import { LSPServerRegistryEntry } from '../../../src/lsp/types';

// ---------------------------------------------------------------------------
// Registry completeness
// ---------------------------------------------------------------------------

describe('LSP Registry', () => {
  const allServers = listBundledLspServers();

  it('contains all local LSP repositories from the OpenQC workspace', () => {
    expect(getBundledLspServerCount()).toBe(17);
    expect(allServers.map(server => server.languageId)).toEqual([
      'abacus',
      'abinit',
      'cif',
      'cp2k',
      'vasp',
      'gaussian',
      'orca',
      'qe',
      'gamess',
      'nwchem',
      'gpumd',
      'gromacs',
      'lammps',
      'mlip',
      'pyatb',
      'pyscf',
      'dpgen',
    ]);
    expect(allServers.map(server => server.id).sort()).toEqual([
      'abacus-lsp',
      'abinit-lsp',
      'cif-lsp',
      'cp2k-lsp-enhanced',
      'dpgen-lsp',
      'gamess-lsp',
      'gaussian-lsp',
      'gpumd-lsp',
      'gromacs-lsp',
      'lammps-lsp',
      'mlip-lsp',
      'nwchem-lsp',
      'orca-lsp',
      'pyatb-lsp',
      'pyscf-lsp',
      'qe-lsp',
      'vasp-lsp',
    ]);
  });

  it('covers every language ID contributed in package.json', () => {
    const languages = packageJson.contributes.languages as Array<{ id: string }>;
    const languageIds = languages.map(l => l.id);

    for (const languageId of languageIds) {
      expect(getLspServerByLanguageId(languageId)).toBeDefined();
    }

    // Registry should not contain language IDs absent from package.json
    const registryLanguageIds = allServers.map(e => e.languageId);
    for (const id of registryLanguageIds) {
      expect(languageIds).toContain(id);
    }
  });

  // ---------------------------------------------------------------------------
  // Stability flags
  // ---------------------------------------------------------------------------

  it('marks the broader local workspace LSPs as experimental', () => {
    const experimentalIds = [
      'abacus',
      'abinit',
      'cif',
      'cp2k',
      'gpumd',
      'gromacs',
      'lammps',
      'mlip',
      'nwchem',
      'pyatb',
      'pyscf',
      'dpgen',
    ];

    for (const languageId of experimentalIds) {
      const entry = getLspServerByLanguageId(languageId);
      expect(entry?.stability).toBe('experimental');
      expect(entry?.defaultBranch).toBeTruthy();
    }
  });

  it('marks the remaining servers as stable', () => {
    const stableIds = ['gaussian', 'orca', 'gamess', 'qe', 'vasp'];
    for (const languageId of stableIds) {
      const entry = getLspServerByLanguageId(languageId);
      expect(entry?.stability).toBe('stable');
      expect(entry?.defaultBranch).toBe('main');
    }
  });

  it('tracks the upstream default branch used for latest-version checks', () => {
    const branches = Object.fromEntries(
      allServers.map(server => [server.id, server.defaultBranch])
    );

    expect(branches).toMatchObject({
      'cif-lsp': 'master',
      'cp2k-lsp-enhanced': 'develop',
      'lammps-lsp': 'master',
    });

    for (const entry of allServers) {
      expect(entry.defaultBranch).toMatch(/^(main|master|develop)$/);
    }
  });

  // ---------------------------------------------------------------------------
  // Lookup helpers
  // ---------------------------------------------------------------------------

  describe('getLspServerByLanguageId', () => {
    it('returns the correct entry for a known language ID', () => {
      const entry = getLspServerByLanguageId('gaussian');
      expect(entry).toMatchObject({
        id: 'gaussian-lsp',
        name: 'Gaussian',
        executable: 'gaussian-lsp',
      });
    });

    it('returns undefined for an unknown language ID', () => {
      expect(getLspServerByLanguageId('nonexistent')).toBeUndefined();
    });
  });

  describe('getLspServerBySoftwareName', () => {
    it('matches human-readable name case-insensitively', () => {
      expect(getLspServerBySoftwareName('Gaussian')).toBeDefined();
      expect(getLspServerBySoftwareName('gaussian')).toBeDefined();
      expect(getLspServerBySoftwareName('GAUSSIAN')).toBeDefined();
    });

    it('matches registry id case-insensitively', () => {
      expect(getLspServerBySoftwareName('gaussian-lsp')).toBeDefined();
      expect(getLspServerBySoftwareName('Gaussian-LSP')).toBeDefined();
    });

    it('returns undefined for an unknown name', () => {
      expect(getLspServerBySoftwareName('unknown')).toBeUndefined();
    });
  });

  // ---------------------------------------------------------------------------
  // Alignment with package.json configuration defaults
  // ---------------------------------------------------------------------------

  it('keeps registry executables aligned with package.json setting defaults', () => {
    const properties = packageJson.contributes.configuration.properties as Record<
      string,
      { default: unknown }
    >;

    for (const entry of allServers) {
      expect(properties[`openqc.lsp.${entry.languageId}.enabled`]?.default).toBe(entry.enabled);
      expect(properties[`openqc.lsp.${entry.languageId}.path`]?.default).toBe(entry.executable);
      expect(properties[`openqc.lsp.${entry.languageId}.command`]?.default).toBe(entry.executable);
      expect(properties[`openqc.lsp.${entry.languageId}.args`]?.default).toEqual(
        entry.args ?? ['--stdio']
      );
    }
  });

  it('keeps registry file patterns aligned with package.json language contributions', () => {
    const languages = packageJson.contributes.languages as Array<{
      id: string;
      extensions?: string[];
      filenames?: string[];
    }>;

    for (const entry of allServers) {
      const language = languages.find(l => l.id === entry.languageId);
      expect(language).toBeDefined();
      expect(language?.extensions || []).toEqual(entry.fileExtensions.map(ext => `.${ext}`));
      expect(language?.filenames || []).toEqual([...entry.fileNames]);
    }
  });

  // ---------------------------------------------------------------------------
  // listBundledLspServers immutability
  // ---------------------------------------------------------------------------

  it('returns a new array on each call (mutations do not affect the registry)', () => {
    const first = listBundledLspServers();
    const second = listBundledLspServers();

    expect(first).not.toBe(second);
    expect(first).toEqual(second);

    first.push({} as LSPServerRegistryEntry);
    expect(listBundledLspServers()).toHaveLength(17);
  });

  it('defines a local launch strategy for every bundled local LSP', () => {
    for (const entry of allServers) {
      expect(entry.localLaunch).toBeDefined();
      expect(entry.localLaunch?.repoName).toBeTruthy();
    }
  });

  it('defines required registry metadata for every bundled local LSP', () => {
    for (const entry of allServers) {
      expect(entry.executable).toBeTruthy();
      expect(entry.languageId).toBeTruthy();
      expect(Array.isArray(entry.fileExtensions)).toBe(true);
      expect(Array.isArray(entry.fileNames)).toBe(true);
      expect(entry.stability).toMatch(/^(stable|experimental)$/);
      expect(entry.repository).toBeTruthy();
      expect(entry.repositoryUrl).toBeTruthy();
      expect(entry.defaultBranch).toBeTruthy();
      expect(entry.localLaunch).toBeDefined();
    }
  });

  it('defines canonical agent CLI operation metadata for every bundled local LSP', () => {
    const expectedOperations = ['check', 'context', 'complete', 'hover', 'symbols', 'fix'];

    for (const entry of allServers) {
      const readiness = getLspDiagnosticReadiness(entry.id);
      expect(readiness?.agentCli).toMatch(/-lsp-tool$/);
      expect(readiness?.agentOperations).toEqual(expectedOperations);
      expect(readiness?.agentCliSmokeStatus).toMatch(/^(pending|passing|failing|unavailable)$/);
      expect(readiness?.richDiagnostics).toBe(true);
    }
  });

  it('keeps docs/LSP_COMPATIBILITY.md aligned with every bundled registry entry', () => {
    const doc = fs.readFileSync(path.join(__dirname, '../../../docs/LSP_COMPATIBILITY.md'), 'utf8');

    for (const entry of allServers) {
      expect(doc).toContain(entry.repository);
      expect(doc).toContain(`\`${entry.languageId}\``);
      expect(doc).toContain(`\`${entry.defaultBranch}\``);
    }

    expect(doc).not.toContain('OpenQuantumChemistry/');
  });

  it('keeps special local launch defaults aligned with sibling repos', () => {
    expect(getLspServerByLanguageId('cp2k')).toMatchObject({
      args: [],
      localLaunch: {
        kind: 'pythonFunction',
        repoName: 'cp2k-lsp-enhanced',
        sourcePath: '.',
      },
    });

    expect(getLspServerByLanguageId('cif')).toMatchObject({
      localLaunch: {
        kind: 'nodeScript',
        repoName: 'cif-lsp',
        scriptPath: 'server/out/server.js',
      },
    });

    expect(getLspServerByLanguageId('lammps')).toMatchObject({
      executable: 'lmp-lsp',
      args: [],
      localLaunch: {
        kind: 'cargoBinary',
        repoName: 'lammps-lsp',
        binaryName: 'lmp-lsp',
      },
    });
  });

  // ---------------------------------------------------------------------------
  // Repository references use newtontech, not OpenQuantumChemistry
  // ---------------------------------------------------------------------------

  it('uses newtontech org in all repository references', () => {
    for (const entry of allServers) {
      expect(entry.repository).toMatch(/^newtontech\//);
      expect(entry.repositoryUrl).toContain('github.com/newtontech/');
    }
  });
});
