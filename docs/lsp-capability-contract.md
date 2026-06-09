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
