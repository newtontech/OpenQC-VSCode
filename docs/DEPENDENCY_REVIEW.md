# Dependency and License Review for Scientific Integrations

## Purpose

This document tracks dependency and license compatibility for the OpenQC VS Code extension's scientific integrations.

## Review Checklist

For each new dependency (npm or Python):

- [ ] **License compatible**: Dependency license is MIT, Apache-2.0, BSD-2/3-Clause, or ISC.
- [ ] **No copyleft**: No GPL/AGPL dependencies in extension runtime (dev-only is acceptable).
- [ ] **Package size**: Evaluate impact on VSIX package size.
- [ ] **Bundle impact**: Dependency does not accidentally ship large native binaries.
- [ ] **Pinned version**: Version is pinned in package.json / requirements.txt.
- [ ] **Security audit**: No known high/critical vulnerabilities in current version.
- [ ] **Maintenance**: Package has recent commits and responds to issues.

## Current npm Dependencies

| Package | Version | License | Purpose | Notes |
|---------|---------|---------|---------|-------|
| 3dmol | ^2.5.5 | MIT | 3D molecular visualization | Bundled in webview |
| ngl | ^2.2.1 | MIT | Protein/molecule visualization | Legacy viewer |
| plotly.js-dist-min | ^3.6.0 | MIT | Data plotting | |
| vscode-languageclient | ^9.0.1 | MIT | LSP client | |

## Optional Python Dependencies

These are NOT bundled with the extension. They are checked at runtime and missing packages produce actionable setup hints.

| Package | License | Purpose | Install Hint |
|---------|---------|---------|-------------|
| pymatgen | MIT | Structure parsing, CIF, POSCAR, supercell | `pip install pymatgen` |
| ase | LGPL-2.1+ | Structure I/O, format conversion | `pip install ase` |
| cclib | MIT | Calculation output parsing | `pip install cclib` |
| dpdata | LGPL-3.0 | Trajectory datasets, ML potentials | `pip install dpdata` |
| spglib | BSD-3-Clause | Space group analysis | `pip install spglib` |
| numpy | BSD-3-Clause | Numerical arrays (transitive) | `pip install numpy` |

### License Notes

- **ASE** (LGPL-2.1+): Used as optional backend only. Extension works without it. LGPL is compatible with MIT extension code when used as a separate process via subprocess bridge.
- **dpdata** (LGPL-3.0): Same as ASE — optional, separate process, no linking.

## External Tools (Optional, Not Bundled)

| Tool | License | Purpose |
|------|---------|---------|
| Multiwfn | Free academic | Wavefunction analysis |
| Open Babel | GPL-2.0 | Format conversion (external process only) |
| c2x | GPL-2.0 | CASTEP density conversion (external process only) |

### Safety Notes

- External analyzers are **disabled by default**.
- User must configure path and confirm before execution.
- GPL tools run as separate processes — no linking with extension code.

## Review Process

1. Add new dependency to the appropriate table above.
2. Run `npm audit` for npm packages.
3. Check license compatibility using SPDX identifiers.
4. Document any exceptions with justification.
5. Update this file in the same PR that adds the dependency.
