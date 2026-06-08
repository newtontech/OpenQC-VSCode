import packageJson from '../../../package.json';
import {
  getBundledLspServerCount,
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

  it('contains exactly seven entries matching the issue spec', () => {
    expect(getBundledLspServerCount()).toBe(7);
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

  it('marks nwchem-lsp and cp2k-lsp-enhanced as experimental', () => {
    const nwchem = getLspServerByLanguageId('nwchem');
    const cp2k = getLspServerByLanguageId('cp2k');

    expect(nwchem?.stability).toBe('experimental');
    expect(nwchem?.defaultBranch).toBe('feature/nwchem-parser');

    expect(cp2k?.stability).toBe('experimental');
    expect(cp2k?.defaultBranch).toBe('develop');
  });

  it('marks the remaining servers as stable', () => {
    const stableIds = ['gaussian', 'orca', 'gamess', 'qe', 'vasp'];
    for (const languageId of stableIds) {
      const entry = getLspServerByLanguageId(languageId);
      expect(entry?.stability).toBe('stable');
      expect(entry?.defaultBranch).toBeUndefined();
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
      expect(language?.extensions || []).toEqual(
        entry.fileExtensions.map(ext => `.${ext}`)
      );
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
    expect(listBundledLspServers()).toHaveLength(7);
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
