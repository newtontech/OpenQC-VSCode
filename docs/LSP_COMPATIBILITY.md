# LSP Compatibility Matrix

This document provides a detailed feature-by-feature status for each standalone language server in the newtontech computational chemistry LSP family, plus its integration status in OpenQC.

Last updated: 2026-06-08 (verified against upstream READMEs and `pyproject.toml` versions)

## Quick Reference

| LSP Server | Version | Domain | File Types | OpenQC Integration |
|------------|---------|--------|------------|--------------------|
| [orca-lsp](https://github.com/newtontech/orca-lsp) | 0.5.4 | ORCA | `.inp` | Syntax + Visualization |
| [gamess-lsp](https://github.com/newtontech/gamess-lsp) | 0.1.0 | GAMESS (US) | `.inp` | Syntax + Visualization |
| [gaussian-lsp](https://github.com/newtontech/gaussian-lsp) | 0.2.11 | Gaussian | `.gjf`, `.com` | Syntax + Visualization |
| [nwchem-lsp](https://github.com/newtontech/nwchem-lsp) | 0.5.0 | NWChem | `.nw`, `.nwinp` | Syntax + Visualization |
| [qe-lsp](https://github.com/newtontech/qe-lsp) | 0.1.0 | Quantum ESPRESSO | `.in`, `.pw.in`, `.relax.in`, `.vc-relax.in`, `.scf.in`, `.nscf.in`, `.bands.in`, `.ph.in`, `.dos.in` | Syntax + Visualization |
| [VASP-LSP](https://github.com/newtontech/VASP-LSP) | 0.4.4 | VASP | `INCAR`, `POSCAR`, `KPOINTS` | Syntax + Visualization |
| [cp2k-lsp-enhanced](https://github.com/newtontech/cp2k-lsp-enhanced) | 0.9.1 | CP2K | `.inp` | Syntax + Visualization |

## Feature Matrix

### Core LSP Capabilities

| Feature | orca-lsp | gamess-lsp | gaussian-lsp | nwchem-lsp | qe-lsp | VASP-LSP | cp2k-lsp-enhanced |
|---------|----------|------------|--------------|------------|--------|----------|--------------------|
| Parser | Full | Full | Full | Full | Full | Full | Full |
| Diagnostics | Yes | Yes | Yes | Yes | Yes | Yes | Yes (cp2klint) |
| Completion | Yes | Yes | Yes | Yes | Yes | Yes | Yes |
| Hover Docs | Yes | Yes | Yes | Yes | Yes | Yes | Yes |
| Formatting | No | Yes | Yes | Yes | No | Yes | No |
| Code Actions | Yes (quick fixes) | Yes (quick fixes) | No | Yes (quick fixes) | No | Yes (quick fixes) | No |
| Document Symbols | No | Yes | No | Yes | No | No | No |
| Folding Ranges | No | No | No | Yes | No | No | No |
| Go to Definition | No | Yes | No | Yes | No | No | No |
| Find References | No | Yes | No | Yes | No | No | No |
| Rename | No | Yes | No | Yes | No | No | No |
| Workspace Symbols | No | Yes | No | Yes | No | No | No |
| Semantic Highlighting | No | No | No | Yes | No | No | No |
| Inlay Hints | No | No | No | Yes | No | No | No |
| Snippet Completions | No | Yes | No | No | No | No | No |

### Status Key

- **Full**: Complete, tested parser for the target format.
- **Yes**: Feature is implemented and active.
- **No**: Feature is not yet implemented.

## Per-Server Details

### orca-lsp (v0.5.4)

- **Repository**: [newtontech/orca-lsp](https://github.com/newtontech/orca-lsp)
- **Description**: LSP implementation for ORCA quantum chemistry software
- **Supported file types**: `.inp` (ORCA input)
- **Parser status**: Full ORCA input file parsing with method/basis/job-type recognition
- **Diagnostics**: Invalid keyword detection, parameter validation, missing required sections, memory and parallelization warnings
- **Completion**: Methods (DFT functionals, wavefunction methods), basis sets (Pople, Karlsruhe, Dunning), job types (SP, OPT, FREQ, etc.), %blocks (%maxcore, %pal, %method, etc.)
- **Hover**: Context-aware documentation for keywords
- **Formatting**: Not implemented
- **Code actions**: Quick fixes for common errors
- **Test command**: `pytest`
- **Build/install**: `pip install orca-lsp` / `pip install -e ".[dev]"`
- **Coverage**: 100% (320 tests)
- **OpenQC integration**: Syntax highlighting, visualization entry points, molecule preview

### gamess-lsp (v0.1.0)

- **Repository**: [newtontech/gamess-lsp](https://github.com/newtontech/gamess-lsp)
- **Description**: LSP implementation for GAMESS (US) input files
- **Supported file types**: `.inp` (GAMESS input)
- **Parser status**: Full parser for GAMESS $group/$END blocks
- **Diagnostics**: Unknown group detection, unclosed section warnings, syntax issues
- **Completion**: $ groups, keywords, value suggestions after `=`
- **Hover**: Inline documentation for GAMESS keywords and groups
- **Formatting**: Automatic formatting with consistent indentation
- **Code actions**: Add missing $END, suggest corrections for unknown groups, add required keywords (e.g., RUNTYP for $CONTRL)
- **Document symbols**: Navigation for $ groups and keywords
- **Go to definition**: Navigate to group or keyword definitions
- **Find references**: Find all occurrences of groups and keywords
- **Rename**: Rename groups and keywords across the document
- **Workspace symbols**: Search for symbols across all open GAMESS files
- **Snippets**: Templates for common GAMESS calculations (water molecule, DFT optimization, MP2, frequency, TD-DFT, etc.)
- **Test command**: `pytest`
- **Build/install**: `pip install gamess-lsp` / `pip install -e ".[dev]"`
- **OpenQC integration**: Syntax highlighting, snippet support, visualization entry points

### gaussian-lsp (v0.2.11)

- **Repository**: [newtontech/gaussian-lsp](https://github.com/newtontech/gaussian-lsp)
- **Description**: LSP implementation for Gaussian quantum chemistry software
- **Supported file types**: `.gjf`, `.com` (Gaussian input)
- **Parser status**: Full parser with route line parsing, periodic table support (118 elements), ModRedundant and ONIOM layer support
- **Diagnostics**: Error and warning detection for Gaussian keywords
- **Completion**: Methods (HF, DFT, post-HF, semi-empirical), basis sets (Pople, Dunning), job types
- **Hover**: Documentation for Gaussian keywords
- **Formatting**: Consistent file structure formatting
- **Code actions**: Not implemented
- **Test command**: `python -m pytest` (Python) / `npm run test:ts` (TypeScript parser)
- **Build/install**: `pip install gaussian-lsp` / `pip install -e ".[dev]"`
- **Additional**: Includes a TypeScript parser under `src/parsers` covered by `npm run typecheck` and `npm run test:ts:coverage`
- **OpenQC integration**: Syntax highlighting, visualization entry points, molecule preview

### nwchem-lsp (v0.5.0)

- **Repository**: [newtontech/nwchem-lsp](https://github.com/newtontech/nwchem-lsp)
- **Description**: LSP implementation for NWChem quantum chemistry software
- **Supported file types**: `.nw`, `.nwinp` (NWChem input)
- **Parser status**: Full parser with section-block recognition
- **Diagnostics**: Unclosed section blocks, unknown basis sets and functionals, invalid task operations, missing required blocks
- **Completion**: Top-level keywords (geometry, basis, scf, dft, etc.), section-specific keywords, basis sets and DFT functionals, task operations and theories
- **Hover**: Inline help for all keywords
- **Formatting**: Automatic code formatting
- **Code actions**: Add missing `end` keywords, remove unexpected `end` keywords, correct common typos (gemoetry -> geometry), add missing `start` directive
- **Document symbols**: Outline view and navigation
- **Folding ranges**: Code folding for sections (v0.5.0)
- **Go to definition**: Navigate from `end` to section start
- **Find references**: Navigate to section occurrences (v0.5.0)
- **Rename**: Rename sections safely (v0.5.0)
- **Workspace symbols**: Global symbol search across all open documents (v0.4.0)
- **Semantic highlighting**: Token-based syntax coloring (v0.4.0)
- **Inlay hints**: Inline information display for units, charge states (v0.4.0)
- **Test command**: `pytest`
- **Build/install**: `pip install nwchem-lsp` / `pip install -e ".[dev]"`
- **Default branch**: `feature/nwchem-parser`
- **Latest release**: v0.2.0 (2026-03-03) — note: pyproject.toml reports 0.5.0 on the default branch
- **OpenQC integration**: Syntax highlighting, visualization entry points, molecule preview

### qe-lsp (v0.1.0)

- **Repository**: [newtontech/qe-lsp](https://github.com/newtontech/qe-lsp)
- **Description**: LSP implementation for Quantum ESPRESSO
- **Supported file types**: `.in`, `.pw.in`, `.relax.in`, `.vc-relax.in`, `.scf.in`, `.nscf.in`, `.bands.in`, `.ph.in`, `.dos.in`
- **Parser status**: Full namelist parsing and validation, card validation (ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS), pseudopotential and element validation, lattice parameter checking
- **Diagnostics**: Error and warning detection for QE namelists and cards
- **Completion**: QE namelists, keywords, and cards
- **Hover**: Documentation for QE keywords and namelists
- **Formatting**: Not implemented
- **Code actions**: Not implemented
- **Test command**: `python -m pytest`
- **Build/install**: `pip install qe-lsp` / `pip install -e ".[dev]"`
- **OpenQC integration**: Syntax highlighting, visualization entry points

### VASP-LSP (v0.4.4)

- **Repository**: [newtontech/VASP-LSP](https://github.com/newtontech/VASP-LSP)
- **Description**: LSP implementation for VASP input/output files
- **Supported file types**: `INCAR`, `POSCAR`, `KPOINTS`
- **Parser status**: Full parser for all three primary VASP input files
- **Diagnostics**: Unknown parameter detection, value type checking, range validation, parameter dependency checks, common configuration warnings
- **Completion**: INCAR parameter names, parameter values (enums, booleans), context-aware suggestions
- **Hover**: Parameter description, valid values/range, default value, related parameters
- **Formatting**: INCAR (parameters grouped by category, aligned values), POSCAR (consistent coordinate precision), KPOINTS (normalized grid types)
- **Code actions**: Add missing SIGMA when ISMEAR >= 0, add missing MAGMOM when ISPIN = 2, add missing LDAU parameters, remove conflicting NPAR/NCORE, fix common tag typos
- **Test command**: `pytest --cov=src/vasp_lsp --cov-report=term-missing`
- **Coverage**: 100%
- **Build/install**: `pip install vasp-lsp` / `pip install -e ".[dev]"`
- **Additional**: Includes VSCode extension in `editors/vscode/`; supports TCP mode for debugging (`vasp-lsp --tcp`); full type annotations
- **OpenQC integration**: Syntax highlighting, visualization entry points, structure preview (POSCAR)

### cp2k-lsp-enhanced (v0.9.1)

- **Repository**: [newtontech/cp2k-lsp-enhanced](https://github.com/newtontech/cp2k-lsp-enhanced)
- **Description**: Fully validating pure-python CP2K input file tools including preprocessing capabilities
- **Supported file types**: `.inp` (CP2K input)
- **Parser status**: Full parser with preprocessing, variable interpolation, @include support, XML-based format specification validation
- **Diagnostics**: Full validation via `cp2klint` and `cp2k-datafile-lint`; XML-backed syntax checking
- **Completion**: Keyword and section completion driven by XML specification
- **Hover**: Contextual help derived from XML specification
- **Formatting**: Not implemented (conversion to/from JSON/YAML preserves structure)
- **Code actions**: Not implemented
- **Additional tools**: `cp2klint` (linter), `fromcp2k` (to JSON/YAML/AiiDA), `tocp2k` (from JSON/YAML), `cp2kgen` (generate variants), `cp2kget` (extract values from restart files)
- **Test command**: `pytest`
- **Server command**: `cp2k-language-server` (via stdio)
- **Build/install**: `pip install cp2k-input-tools` (core) / `pip install cp2k-input-tools[lsp]` (with LSP) / `pip install cp2k-input-tools[yaml]` (with YAML)
- **Default branch**: `develop`
- **OpenQC integration**: Syntax highlighting, visualization entry points, CP2K parser/linter tooling alignment

## Build and Test Commands

### Standalone LSP Servers

All servers use the pygls framework and share a common build/test pattern:

| Step | Command |
|------|---------|
| Install (user) | `pip install <package-name>` |
| Install (dev) | `pip install -e ".[dev]"` |
| Run server | `<server-name>` (stdio) |
| Run tests | `pytest` |
| Run with coverage | `pytest --cov --cov-report=term-missing` |
| Format check | `black --check src/ tests/` |
| Lint | `ruff check src/ tests/` or `flake8 src/ tests/` |
| Type check | `mypy src/` |

#### Server commands

| Server | `pip install` name | Server command |
|--------|-------------------|----------------|
| orca-lsp | `orca-lsp` | `orca-lsp` |
| gamess-lsp | `gamess-lsp` | `gamess-lsp` |
| gaussian-lsp | `gaussian-lsp` | `gaussian-lsp` |
| nwchem-lsp | `nwchem-lsp` | `nwchem-lsp` |
| qe-lsp | `qe-lsp` | `qe-lsp` |
| VASP-LSP | `vasp-lsp` | `vasp-lsp --stdio` (also `--tcp` for debugging) |
| cp2k-lsp-enhanced | `cp2k-input-tools[lsp]` | `cp2k-language-server` |

### OpenQC-VSCode

| Step | Command |
|------|---------|
| Install | `npm install` |
| Compile | `npm run compile` |
| Lint | `npm run lint` |
| Type check | `npm run typecheck` |
| Format check | `npm run format:check` |
| Unit tests | `npm run test:unit` |
| Integration tests | `npm run test:integration` |
| Full CI (PR gate) | `npm run ci:pr` |
| Package VSIX | `npm run package:vsix` |

## VSCode / OpenQC Integration Status

| LSP Server | Language ID | Syntax Grammar | Snippets | Visualization | Sidebar Integration | Direct LSP Wiring |
|------------|-------------|----------------|----------|---------------|---------------------|-------------------|
| orca-lsp | `orca` | Yes | No | Yes (3D molecule) | Yes (Molecules) | OpenQC modules |
| gamess-lsp | `gamess` | Yes | Yes | Yes (3D molecule) | Yes (Molecules) | OpenQC modules |
| gaussian-lsp | `gaussian` | Yes | No | Yes (3D molecule) | Yes (Molecules) | OpenQC modules |
| nwchem-lsp | `nwchem` | Yes | No | Yes (3D molecule) | Yes (Molecules) | OpenQC modules |
| qe-lsp | `quantum-espresso` | Yes | No | Yes (3D structure) | Yes (Molecules) | OpenQC modules |
| VASP-LSP | `vasp` (INCAR/POSCAR/KPOINTS) | Yes | No | Yes (3D structure) | Yes (Molecules) | OpenQC modules |
| cp2k-lsp-enhanced | `cp2k` | Yes | No | Yes (3D molecule) | Yes (Molecules) | OpenQC modules |

### Integration Notes

- **Direct LSP Wiring**: All servers currently integrate through OpenQC's internal modules rather than spawning standalone LSP processes. The long-term goal is to support direct `stdio` connections to each standalone server.
- **Syntax Grammar**: All formats have TextMate grammars in OpenQC's `syntaxes/` directory.
- **Visualization**: Formats that carry molecular or crystal structure information (coordinates) support the 3D viewer via the editor toolbar icon.
- **Sidebar**: The Molecules sidebar tracks parsed structures from supported file types.

## Editor Integration (Standalone)

All standalone servers communicate via stdio and can be used with any LSP-capable editor. Below are tested configurations:

### Neovim (nvim-lspconfig)

```lua
-- Example for gamess-lsp
local lspconfig = require('lspconfig')
lspconfig.gamess.setup {
  cmd = {"gamess-lsp"},
  filetypes = {"gamess"},
  root_dir = lspconfig.util.root_pattern("*.inp"),
}

-- Example for nwchem-lsp
lspconfig.nwchem.setup {
  cmd = {"nwchem-lsp"},
  filetypes = {"nw", "nwinp"},
}

-- Example for VASP-LSP
require'lspconfig'.vasp_lsp.setup{}
```

### Emacs (lsp-mode)

```elisp
;; Example for gamess-lsp
(lsp-register-client
 (make-lsp-client :new-connection (lsp-stdio-connection "gamess-lsp")
                  :major-modes '(gamess-mode)
                  :server-id 'gamess-lsp))

;; Example for nwchem-lsp
(lsp-register-client
 (make-lsp-client :new-connection (lsp-stdio-connection "nwchem-lsp")
                  :major-modes '(nwchem-mode)
                  :server-id 'nwchem-lsp))
```

### VS Code (without OpenQC)

Add to `settings.json`:

```json
{
  "languageserver": {
    "gamess": {
      "command": "gamess-lsp",
      "filetypes": ["gamess"],
      "rootPatterns": ["*.inp"]
    },
    "nwchem": {
      "command": "nwchem-lsp",
      "filetypes": ["nw", "nwinp"],
      "rootPatterns": ["*.nw", "*.nwinp"]
    }
  }
}
```

VASP-LSP ships its own VS Code extension under `editors/vscode/`.

## Alignment Rules

1. Keep file extension and language-id decisions consistent across OpenQC and each standalone LSP.
2. Prefer shared fixtures for parser behavior whenever the format appears in both places.
3. Document any OpenQC-only behavior, such as visualization commands or sidebar integration.
4. Track mismatches with the `lsp-alignment` label.
5. Release notes should state which standalone LSP versions or commit SHAs were smoke tested.

## Smoke Test Checklist

Before a Marketplace release, open one minimal sample for each format and verify:

- [ ] Syntax highlighting activates.
- [ ] The expected language server or parser path starts.
- [ ] Diagnostics appear for one intentionally invalid input.
- [ ] Hover or completion works for one common keyword.
- [ ] Visualization entry points do not break for molecule or structure-bearing inputs.
