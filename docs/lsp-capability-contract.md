# LSP Capability Contract

This document defines the capability contract and compatibility matrix for all OpenQC language servers.

## Supported Language Servers

| Server | Language | Package | Version |
|--------|----------|---------|---------|
| nwchem-lsp | NWChem | `nwchem-lsp` | latest |
| qe-lsp | Quantum ESPRESSO | `qe-lsp` | latest |
| gaussian-lsp | Gaussian | `gaussian-lsp` | latest |
| orca-lsp | ORCA | `orca-lsp` | latest |
| gamess-lsp | GAMESS | `gamess-lsp` | latest |
| cp2k-lsp-enhanced | CP2K | `cp2k-lsp-enhanced` | latest |
| abacus-lsp | ABACUS | `abacus-lsp` | experimental |
| abinit-lsp | ABINIT | `abinit-lsp` | experimental |
| gpumd-lsp | GPUMD | `gpumd-lsp` | experimental |
| gromacs-lsp | GROMACS | `gromacs-lsp` | experimental |
| lammps-lsp | LAMMPS | `lammps-lsp` (cargo) | experimental |
| pyscf-lsp | PySCF | `pyscf-lsp` | experimental |
| pyatb-lsp | PyATB | `pyatb-lsp` | experimental |
| mlip-lsp | MLIP | `mlip-lsp` | experimental |
| cif-lsp | CIF | `cif-lsp` (npm) | experimental |

## LSP Feature Matrix

| Feature | LSP Method | nwchem | qe | gaussian | orca | gamess |
|---------|-----------|--------|----|----------|------|--------|
| Diagnostics | textDocument/diagnostic | ✅ | ✅ | ✅ | ✅ | ✅ |
| Completion | textDocument/completion | ✅ | ✅ | ✅ | ✅ | ✅ |
| Hover | textDocument/hover | ✅ | ✅ | ✅ | ✅ | ✅ |
| Definition | textDocument/definition | ✅ | ✅ | ✅ | ✅ | ✅ |
| References | textDocument/references | ✅ | ✅ | ✅ | ✅ | ✅ |
| Formatting | textDocument/formatting | ✅ | ✅ | ✅ | ✅ | ✅ |
| Range Formatting | textDocument/rangeFormatting | ✅ | ✅ | ✅ | ✅ | ✅ |
| Code Actions | textDocument/codeAction | ✅ | ✅ | ✅ | ✅ | ✅ |
| Rename | textDocument/rename | ✅ | ✅ | ✅ | ✅ | ✅ |
| Document Symbols | textDocument/documentSymbol | ✅ | ✅ | ✅ | ✅ | ✅ |
| Lint Checks | Custom diagnostics | ✅ | ✅ | ✅ | ✅ | ✅ |
| Type Checking | Custom diagnostics | ✅ | ✅ | ✅ | ✅ | ✅ |

## Custom Commands

| Command | Description |
|---------|-------------|
| `{lang}.diagnosticsSnapshot` | Returns JSON snapshot of current diagnostics |

## Client Requirements

- VS Code >= 1.85
- Python >= 3.10
- pygls >= 2.0

## Activation Events

All LSP servers activate on file extensions:
- `.nw`, `.nwi` → nwchem-lsp
- `.in`, `.qe` → qe-lsp
- `.gjf`, `.com` → gaussian-lsp
- `.inp` → orca-lsp, gamess-lsp
- `INPUT`, `STRU`, `KPT` → abacus-lsp
- `.abi`, `.abinit` → abinit-lsp
- `run.in`, `nep.in` → gpumd-lsp
- `.top`, `.itp`, `.mdp`, `.gro` → gromacs-lsp
- `.lmp`, `.lammps`, `.lmps`, `in.lammps` → lammps-lsp
- `.pyscf.py`, `run_pyscf.py` → pyscf-lsp
- `.pyatb.py`, `run_pyatb.py` → pyatb-lsp
- `.mlip.json`, `.mlip.yaml`, `.mlip.yml` → mlip-lsp
- `.cif` → cif-lsp
