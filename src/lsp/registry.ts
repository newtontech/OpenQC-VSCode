/**
 * Bundled LSP Server Registry
 *
 * Static, version-controlled registry that serves as the single source of truth
 * for known LSP servers. Replaces runtime GitHub discovery with a bundled lookup.
 *
 * @module lsp/registry
 * @see https://github.com/newtontech/OpenQC-VSCode/issues/96
 */

import { LSPServerRegistryEntry } from './types';

// ---------------------------------------------------------------------------
// Registry data
// ---------------------------------------------------------------------------

/**
 * The canonical list of bundled LSP server entries.
 *
 * Order does not matter — lookups are by ID or languageId.
 */
const BUNDLED_LSP_SERVERS: readonly LSPServerRegistryEntry[] = [
  {
    id: 'gaussian-lsp',
    name: 'Gaussian',
    repository: 'newtontech/gaussian-lsp',
    executable: 'gaussian-lsp',
    languageId: 'gaussian',
    fileExtensions: ['gjf', 'com'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/gaussian-lsp',
    stability: 'stable',
  },
  {
    id: 'orca-lsp',
    name: 'ORCA',
    repository: 'newtontech/orca-lsp',
    executable: 'orca-lsp',
    languageId: 'orca',
    fileExtensions: ['inp'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/orca-lsp',
    stability: 'stable',
  },
  {
    id: 'gamess-lsp',
    name: 'GAMESS',
    repository: 'newtontech/gamess-lsp',
    executable: 'gamess-lsp',
    languageId: 'gamess',
    fileExtensions: ['inp'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/gamess-lsp',
    stability: 'stable',
  },
  {
    id: 'qe-lsp',
    name: 'Quantum ESPRESSO',
    repository: 'newtontech/qe-lsp',
    executable: 'qe-lsp',
    languageId: 'qe',
    fileExtensions: [
      'in',
      'pw.in',
      'relax.in',
      'vc-relax.in',
      'scf.in',
      'nscf.in',
      'bands.in',
      'ph.in',
      'dos.in',
    ],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/qe-lsp',
    stability: 'stable',
  },
  {
    id: 'vasp-lsp',
    name: 'VASP',
    repository: 'newtontech/VASP-LSP',
    executable: 'vasp-lsp',
    languageId: 'vasp',
    fileExtensions: [],
    fileNames: [
      'INCAR',
      'POSCAR',
      'KPOINTS',
      'POTCAR',
      'CONTCAR',
      'OSZICAR',
      'OUTCAR',
      'vasprun.xml',
    ],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/VASP-LSP',
    stability: 'stable',
  },
  {
    id: 'nwchem-lsp',
    name: 'NWChem',
    repository: 'newtontech/nwchem-lsp',
    executable: 'nwchem-lsp',
    languageId: 'nwchem',
    fileExtensions: ['nw', 'nwinp'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/nwchem-lsp',
    stability: 'experimental',
    defaultBranch: 'feature/nwchem-parser',
  },
  {
    id: 'cp2k-lsp-enhanced',
    name: 'CP2K',
    repository: 'newtontech/cp2k-lsp-enhanced',
    executable: 'cp2k-language-server',
    languageId: 'cp2k',
    fileExtensions: ['inp'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/cp2k-lsp-enhanced',
    stability: 'experimental',
    defaultBranch: 'develop',
  },
] as const;

// ---------------------------------------------------------------------------
// Lookup helpers
// ---------------------------------------------------------------------------

/**
 * Return the registry entry for the given VS Code language ID, or `undefined`.
 */
export function getLspServerByLanguageId(languageId: string): LSPServerRegistryEntry | undefined {
  return BUNDLED_LSP_SERVERS.find(entry => entry.languageId === languageId);
}

/**
 * Return the registry entry for the given software name (case-insensitive),
 * or `undefined`.
 *
 * Matches against both the human-readable `name` ("Gaussian") and the
 * registry `id` ("gaussian-lsp").
 */
export function getLspServerBySoftwareName(name: string): LSPServerRegistryEntry | undefined {
  const lowered = name.toLowerCase();
  return BUNDLED_LSP_SERVERS.find(
    entry => entry.name.toLowerCase() === lowered || entry.id.toLowerCase() === lowered
  );
}

/**
 * Return a shallow copy of the full bundled registry.
 *
 * The returned array is safe to mutate; the internal registry is frozen.
 */
export function listBundledLspServers(): LSPServerRegistryEntry[] {
  return [...BUNDLED_LSP_SERVERS];
}

/**
 * Return the number of entries in the registry (useful for tests).
 */
export function getBundledLspServerCount(): number {
  return BUNDLED_LSP_SERVERS.length;
}
