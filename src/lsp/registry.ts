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
 * Order follows the package.json language contribution order.
 */
const BUNDLED_LSP_SERVERS: readonly LSPServerRegistryEntry[] = [
  {
    id: 'abacus-lsp',
    name: 'ABACUS',
    repository: 'newtontech/abacus-lsp',
    executable: 'abacus-lsp',
    languageId: 'abacus',
    fileExtensions: [],
    fileNames: ['INPUT', 'STRU', 'KPT'],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/abacus-lsp',
    stability: 'experimental',
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'abacus-lsp',
      importPath: 'abacus_lsp.cli',
      functionName: 'lsp_main',
    },
  },
  {
    id: 'abinit-lsp',
    name: 'ABINIT',
    repository: 'newtontech/abinit-lsp',
    executable: 'abinit-lsp',
    languageId: 'abinit',
    fileExtensions: ['abi', 'abinit'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/abinit-lsp',
    stability: 'experimental',
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'abinit-lsp',
      importPath: 'abinit_lsp.cli',
      functionName: 'lsp_main',
    },
  },
  {
    id: 'cif-lsp',
    name: 'CIF',
    repository: 'newtontech/cif-lsp',
    executable: 'cif-lsp',
    languageId: 'cif',
    fileExtensions: ['cif'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/cif-lsp',
    stability: 'experimental',
    defaultBranch: 'master',
    localLaunch: {
      kind: 'nodeScript',
      repoName: 'cif-lsp',
      scriptPath: 'server/out/server.js',
    },
  },
  {
    id: 'cp2k-lsp-enhanced',
    name: 'CP2K',
    repository: 'newtontech/cp2k-lsp-enhanced',
    executable: 'cp2k-language-server',
    args: [],
    languageId: 'cp2k',
    fileExtensions: ['inp'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/cp2k-lsp-enhanced',
    stability: 'experimental',
    defaultBranch: 'develop',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'cp2k-lsp-enhanced',
      importPath: 'cp2k_input_tools.cli.lsp',
      functionName: 'cp2k_language_server',
      sourcePath: '.',
    },
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
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'VASP-LSP',
      importPath: 'vasp_lsp.server',
      functionName: 'main',
    },
  },
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
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'gaussian-lsp',
      importPath: 'gaussian_lsp.server',
      functionName: 'main',
    },
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
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'orca-lsp',
      importPath: 'orca_lsp.server',
      functionName: 'main',
    },
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
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'qe-lsp',
      importPath: 'qe_lsp.server',
      functionName: 'main',
    },
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
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'gamess-lsp',
      importPath: 'gamess_lsp.server',
      functionName: 'main',
    },
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
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'nwchem-lsp',
      importPath: 'nwchem_lsp.server',
      functionName: 'main',
    },
  },
  {
    id: 'gpumd-lsp',
    name: 'GPUMD',
    repository: 'newtontech/gpumd-lsp',
    executable: 'gpumd-lsp',
    languageId: 'gpumd',
    fileExtensions: [],
    fileNames: ['run.in', 'nep.in'],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/gpumd-lsp',
    stability: 'experimental',
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'gpumd-lsp',
      importPath: 'gpumd_lsp.cli',
      functionName: 'lsp_main',
    },
  },
  {
    id: 'gromacs-lsp',
    name: 'GROMACS',
    repository: 'newtontech/gromacs-lsp',
    executable: 'gromacs-lsp',
    languageId: 'gromacs',
    fileExtensions: ['top', 'itp', 'mdp', 'gro'],
    fileNames: [],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/gromacs-lsp',
    stability: 'experimental',
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'gromacs-lsp',
      importPath: 'gromacs_lsp.cli',
      functionName: 'lsp_main',
      sourcePath: '.',
    },
  },
  {
    id: 'lammps-lsp',
    name: 'LAMMPS',
    repository: 'newtontech/lammps-lsp',
    executable: 'lmp-lsp',
    args: [],
    languageId: 'lammps',
    fileExtensions: ['lmp', 'lammps', 'lmps'],
    fileNames: ['in.lammps'],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/lammps-lsp',
    stability: 'experimental',
    defaultBranch: 'master',
    localLaunch: {
      kind: 'cargoBinary',
      repoName: 'lammps-lsp',
      binaryName: 'lmp-lsp',
    },
  },
  {
    id: 'mlip-lsp',
    name: 'MLIP',
    repository: 'newtontech/mlip-lsp',
    executable: 'mlip-lsp',
    languageId: 'mlip',
    fileExtensions: ['mlip.json', 'mlip.yaml', 'mlip.yml'],
    fileNames: ['mlip.json', 'mlip.yaml', 'mlip.yml'],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/mlip-lsp',
    stability: 'experimental',
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'mlip-lsp',
      importPath: 'mlip_lsp.cli',
      functionName: 'lsp_main',
    },
  },
  {
    id: 'pyatb-lsp',
    name: 'PyATB',
    repository: 'newtontech/pyatb-lsp',
    executable: 'pyatb-lsp',
    languageId: 'pyatb',
    fileExtensions: ['pyatb.py'],
    fileNames: ['run_pyatb.py'],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/pyatb-lsp',
    stability: 'experimental',
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'pyatb-lsp',
      importPath: 'pyatb_lsp.cli',
      functionName: 'lsp_main',
    },
  },
  {
    id: 'pyscf-lsp',
    name: 'PySCF',
    repository: 'newtontech/pyscf-lsp',
    executable: 'pyscf-lsp',
    languageId: 'pyscf',
    fileExtensions: ['pyscf.py'],
    fileNames: ['run_pyscf.py'],
    enabled: true,
    repositoryUrl: 'https://github.com/newtontech/pyscf-lsp',
    stability: 'experimental',
    defaultBranch: 'main',
    localLaunch: {
      kind: 'pythonFunction',
      repoName: 'pyscf-lsp',
      importPath: 'pyscf_lsp.cli',
      functionName: 'lsp_main',
    },
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
