# LSP Capability Contract

This document defines the client-side compatibility contract that OpenQC applies to every bundled language server. Per-server capability details belong in each standalone LSP repository; OpenQC only promises consistent discovery, configuration, startup, and editor integration.

## Supported Language Servers

| Server | Language ID | Package/Command | Default Branch |
|--------|-------------|-----------------|----------------|
| abacus-lsp | `abacus` | `abacus-lsp` | `main` |
| abinit-lsp | `abinit` | `abinit-lsp` | `main` |
| cif-lsp | `cif` | `cif-lsp` | `master` |
| cp2k-lsp-enhanced | `cp2k` | `cp2k-language-server` | `develop` |
| VASP-LSP | `vasp` | `vasp-lsp` | `main` |
| gaussian-lsp | `gaussian` | `gaussian-lsp` | `main` |
| orca-lsp | `orca` | `orca-lsp` | `main` |
| qe-lsp | `qe` | `qe-lsp` | `main` |
| gamess-lsp | `gamess` | `gamess-lsp` | `main` |
| nwchem-lsp | `nwchem` | `nwchem-lsp` | `main` |
| gpumd-lsp | `gpumd` | `gpumd-lsp` | `main` |
| gromacs-lsp | `gromacs` | `gromacs-lsp` | `main` |
| lammps-lsp | `lammps` | `lmp-lsp` | `master` |
| mlip-lsp | `mlip` | `mlip-lsp` | `main` |
| pyatb-lsp | `pyatb` | `pyatb-lsp` | `main` |
| pyscf-lsp | `pyscf` | `pyscf-lsp` | `main` |

## Client-Side Guarantees

| Area | Contract |
|------|----------|
| Activation | OpenQC activates for every language ID listed in `src/lsp/registry.ts`. |
| Startup | OpenQC starts each server over stdio using the configured command and args. |
| Configuration | Every server exposes `enabled`, `path`, `command`, `args`, and `env` settings. |
| Workspace lifecycle | OpenQC keeps one client per language per workspace folder and stops it when no matching documents remain. |
| Latest tracking | Every registry entry records the upstream default branch used for latest-version checks and release smoke tests. |
| Compatibility report | `OpenQC LSP: Generate Compatibility Report` reports the configured command, repository, default branch, stability, and current session status. |

## Optional Server Capabilities

OpenQC does not assume that all servers implement the same optional LSP methods. The extension should tolerate missing support for:

- `textDocument/formatting`
- `textDocument/codeAction`
- `textDocument/definition`
- `textDocument/references`
- `textDocument/rename`
- `textDocument/documentSymbol`
- semantic tokens
- inlay hints

When a server adds or removes one of these methods, update `docs/LSP_COMPATIBILITY.md` and release notes after smoke testing the latest server revision.
