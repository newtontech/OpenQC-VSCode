# LSP Alignment

OpenQC is the VS Code product surface for the newtontech computational chemistry LSP family. Standalone LSP repositories remain useful for editor-agnostic language support; OpenQC should integrate them with a consistent VS Code experience.

> **Full feature-by-feature compatibility matrix**: See [LSP_COMPATIBILITY.md](./LSP_COMPATIBILITY.md) for per-server parser status, diagnostics, completion, hover, formatting, code actions, build commands, and integration details.

## Canonical repositories

| Repository | Domain | OpenQC expectation |
|------------|--------|--------------------|
| `newtontech/gaussian-lsp` | Gaussian `.gjf`, `.com` | Match diagnostics, completion vocabulary, and sample fixtures |
| `newtontech/orca-lsp` | ORCA `.inp` | Match diagnostics, completion vocabulary, and sample fixtures |
| `newtontech/gamess-lsp` | GAMESS `.inp` | Match group/keyword diagnostics, snippets, and fixtures |
| `newtontech/qe-lsp` | Quantum ESPRESSO `.in` variants | Match diagnostics, completion vocabulary, and file detection |
| `newtontech/nwchem-lsp` | NWChem `.nw`, `.nwinp` | Match section parsing, diagnostics, and inlay hints |
| `newtontech/VASP-LSP` | VASP `INCAR`, `POSCAR`, `KPOINTS` | Match parameter validation, quick fixes, and formatting |
| `newtontech/cp2k-lsp-enhanced` | CP2K `.inp` | Reuse or align with CP2K parser/linter behavior |

## Alignment rules

- Keep file extension and language-id decisions consistent across OpenQC and each standalone LSP.
- Prefer shared fixtures for parser behavior whenever the format appears in both places.
- Document any OpenQC-only behavior, such as visualization commands or sidebar integration.
- Track mismatches with the `lsp-alignment` label.
- Release notes should state which standalone LSP versions or commit SHAs were smoke tested.

## Smoke test matrix

Before a Marketplace release, open one minimal sample for each format and verify:

- Syntax highlighting activates.
- The expected language server or parser path starts.
- Diagnostics appear for one intentionally invalid input.
- Hover or completion works for one common keyword.
- Visualization entry points do not break for molecule or structure-bearing inputs.
