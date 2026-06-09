# Dependency Baseline Review

## Overview

This document provides a baseline review of OpenQC-VSCode dependencies after the LSP cleanup and modernization effort.

## Runtime Dependencies

| Package | Version | Purpose | License |
|---------|---------|---------|---------|
| vscode-languageclient | ^9.0.1 | LSP client protocol | MIT |
| ngl | ^2.2.1 | 3D molecular visualization | MIT |
| ws | ^8.0.0 | WebSocket client | MIT |

## Dev Dependencies

| Package | Version | Purpose | License |
|---------|---------|---------|---------|
| @types/vscode | ^1.85.0 | VS Code API types | MIT |
| typescript | ^5.3.0 | Type checking | Apache-2.0 |
| jest | ^29.7.0 | Testing | MIT |
| eslint | ^8.56.0 | Linting | MIT |
| prettier | ^3.2.0 | Formatting | MIT |

## LSP Server Dependencies (Python)

Each LSP server uses:
- `pygls` >= 2.0 — Language Server Protocol framework
- `lsprotocol` — LSP type definitions
- `pytest` >= 7.0 — Testing

## Security Review

- No known CVEs in current versions
- All dependencies use permissive licenses (MIT, Apache-2.0)
- `vscode-languageclient` v10 available but not yet tested for compatibility

## Recommendations

1. Test `vscode-languageclient` v10 upgrade
2. Review `ngl` v2.4.0 for viewer improvements
3. Add `npm audit` to CI pipeline
4. Consider adding Dependabot for automated updates
