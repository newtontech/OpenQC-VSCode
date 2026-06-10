# LSP Family Interface

This document describes the code-sharing architecture across the newtontech quantum chemistry LSP family.

## LSP Repositories

| Repository | Language ID | Executable | File Patterns |
|-----------|-------------|-----------|---------------|
| [vasp-lsp](https://github.com/newtontech/VASP-LSP) | `vasp` | `vasp-lsp` | INCAR, POSCAR, KPOINTS, POTCAR, CONTCAR |
| [gaussian-lsp](https://github.com/newtontech/gaussian-lsp) | `gaussian` | `gaussian-lsp` | *.gjf, *.com |
| [orca-lsp](https://github.com/newtontech/orca-lsp) | `orca` | `orca-lsp` | *.inp |
| [cp2k-lsp-enhanced](https://github.com/newtontech/cp2k-lsp-enhanced) | `cp2k` | `cp2k-language-server` | *.inp |
| [qe-lsp](https://github.com/newtontech/qe-lsp) | `qe` | `qe-lsp` | *.in, *.pw.in, *.scf.in, etc. |
| [gamess-lsp](https://github.com/newtontech/gamess-lsp) | `gamess` | `gamess-lsp` | *.inp |
| [nwchem-lsp](https://github.com/newtontech/nwchem-lsp) | `nwchem` | `nwchem-lsp` | *.nw, *.nwinp |
| [abacus-lsp](https://github.com/newtontech/abacus-lsp) | `abacus` | `abacus-lsp` | INPUT, STRU, KPT |
| [abinit-lsp](https://github.com/newtontech/abinit-lsp) | `abinit` | `abinit-lsp` | *.abi, *.abinit |
| [gpumd-lsp](https://github.com/newtontech/gpumd-lsp) | `gpumd` | `gpumd-lsp` | run.in, nep.in |
| [gromacs-lsp](https://github.com/newtontech/gromacs-lsp) | `gromacs` | `gromacs-lsp` | *.top, *.itp, *.mdp, *.gro |
| [lammps-lsp](https://github.com/newtontech/lammps-lsp) | `lammps` | `lmp-lsp` | *.lmp, *.lammps, *.lmps, in.lammps |
| [pyscf-lsp](https://github.com/newtontech/pyscf-lsp) | `pyscf` | `pyscf-lsp` | *.pyscf.py, run_pyscf.py |
| [pyatb-lsp](https://github.com/newtontech/pyatb-lsp) | `pyatb` | `pyatb-lsp` | *.pyatb.py, run_pyatb.py |
| [mlip-lsp](https://github.com/newtontech/mlip-lsp) | `mlip` | `mlip-lsp` | *.mlip.json, *.mlip.yaml, *.mlip.yml |
| [cif-lsp](https://github.com/newtontech/cif-lsp) | `cif` | `cif-lsp` | *.cif |

## Discovery and Configuration

`LSPDiscovery` in `src/utils/LSPDiscovery.ts` maintains a canonical list of LSP server definitions. On startup it:

1. Checks workspace/global state cache (1-hour TTL)
2. Falls back to GitHub API (`/orgs/OpenQuantumChemistry/repos`)
3. Falls back to hardcoded defaults in `DEFAULT_LSP_SERVER_DEFINITIONS`

Users can override the executable path and enable/disable each LSP via VS Code settings:

```
openqc.lsp.<languageId>.enabled  (boolean)
openqc.lsp.<languageId>.path     (string)
```

## Client Lifecycle

`LSPManager` in `src/managers/LSPManager.ts` manages LSP clients:

- **Workspace-scoped**: One client per language per workspace folder
- **Deduplication**: Concurrent starts for the same language+workspace are coalesced
- **Reference counting**: Clients are stopped when the last matching document closes
- **Graceful disposal**: All clients are stopped on extension deactivation

## Shared Parser Infrastructure

The OpenQC-VSCode extension itself provides shared parsing logic that could inform the individual LSPs:

- **`BaseParser`** (`src/parsers/base.ts`): Abstract parser with shared utilities for comment stripping, key-value parsing, and coordinate extraction
- **`CoordinateCapable`** interface: Uniform contract for parsers that extract atomic coordinates
- **`FileTypeDetector`** (`src/managers/FileTypeDetector.ts`): Content-based format detection used to route files to the correct LSP

## Adding a New LSP

To add support for a new quantum chemistry LSP:

1. Add the LSP server definition to `DEFAULT_LSP_SERVER_DEFINITIONS` in `src/utils/LSPDiscovery.ts`
2. Add language configuration in `package.json` under `contributes.languages`
3. Add TextMate grammar in `package.json` under `contributes.grammars`
4. Add LSP settings in `package.json` under `contributes.configuration.properties`
5. Create a parser in `src/parsers/` extending `BaseParser`
6. Register the parser in `src/parsers/index.ts` `createParser` function
7. Add format detection patterns in `FileTypeDetector`

## Cross-LSP Code Sharing Opportunities

The following modules are candidates for extraction into a shared package that all LSPs could consume:

| Module | Current Location | Sharing Potential |
|--------|-----------------|-------------------|
| Coordinate parsing | `BaseParser.parseCoordinateLine` | High - all parsers extract xyz coordinates |
| Comment stripping | `BaseParser.stripComments` | Medium - comment characters vary by format |
| Keyword validation | Per-parser `validate()` | Medium - common validation patterns |
| Element data | `src/visualizers/atomicData.ts` | High - CPK colors, covalent/vdW radii |
| File type detection | `FileTypeDetector` | Medium - content heuristics could be shared |
