# LSP Compatibility Matrix

This document tracks the feature status of all NewtonTech LSP (Language Server Protocol) implementations. OpenQC-VSCode acts as the hub and integrates with all listed LSP servers.

## Feature Matrix

| Feature | ORCA | GAMESS | Gaussian | NWChem | QE | VASP | CP2K |
|---------|------|--------|----------|--------|-----|------|------|
| **Parser** | ✅ Active | ✅ Active | ✅ Active | ✅ Active (Python) | ✅ Active | ✅ Active | ✅ Active |
| **Diagnostics** | ✅ | ✅ Warnings + errors | ✅ Errors + warnings | ✅ Real-time | ✅ Errors + warnings | ✅ Real-time | ✅ Full validation |
| **Completion** | ✅ Auto | ✅ + Snippets | ✅ Keywords | ✅ Context-aware | ✅ Namelists + cards | ✅ Tags + values | ✅ Via cp2k-lsp |
| **Hover** | ✅ Context-aware | ✅ Keywords + groups | ✅ Keywords | ✅ Inline help | ✅ Keywords | ✅ Parameters | ✅ Via cp2k-lsp |
| **Formatting** | ❌ | ✅ Auto-indent | ✅ | ✅ Auto | ❌ | ✅ INCAR/POSCAR/KPOINTS | ❌ |
| **Code Actions** | ✅ Quick fixes | ✅ Quick fixes | ❌ | ✅ Quick fixes | ❌ | ✅ | ❌ |
| **Inlay Hints** | ❌ | ❌ | ❌ | ✅ Units, charge | ❌ | ❌ | ❌ |
| **Semantic Validation** | ✅ | ✅ Mutual exclusivity | ❌ | ✅ Mutual exclusivity | ❌ | ✅ INCAR/POTCAR/POSCAR | ✅ RUN_TYPE/MOTION |

## Integration Status

| LSP Server | VS Code Extension | OpenQC Integration | Language |
|------------|------------------|-------------------|----------|
| `orca-lsp` | Via OpenQC | ✅ Aligned | Python |
| `gamess-lsp` | Via OpenQC | ✅ Aligned | Python |
| `gaussian-lsp` | Via OpenQC | ✅ Aligned | TypeScript + Python |
| `nwchem-lsp` | Via OpenQC | ✅ Aligned | Python |
| `qe-lsp` | Via OpenQC | ✅ Aligned | Python |
| `VASP-LSP` | Via OpenQC | ✅ Aligned | Python |
| `cp2k-lsp-enhanced` | Via OpenQC | ✅ Aligned | Python |
| `abacus-lsp` | Via OpenQC | 🔧 Experimental | Python |
| `abinit-lsp` | Via OpenQC | 🔧 Experimental | Python |
| `gpumd-lsp` | Via OpenQC | 🔧 Experimental | Python |
| `gromacs-lsp` | Via OpenQC | 🔧 Experimental | Python |
| `lammps-lsp` | Via OpenQC | 🔧 Experimental | Rust |
| `pyscf-lsp` | Via OpenQC | 🔧 Experimental | Python |
| `pyatb-lsp` | Via OpenQC | 🔧 Experimental | Python |
| `mlip-lsp` | Via OpenQC | 🔧 Experimental | Python |
| `cif-lsp` | Via OpenQC | 🔧 Experimental | TypeScript |

## Test Commands

| Repo | Test Command | Build Command |
|------|-------------|---------------|
| `orca-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `gamess-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `gaussian-lsp` | `uv run pytest tests/ -v` | `pip install -e .` |
| `nwchem-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `qe-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `VASP-LSP` | `python -m pytest tests/ -v` | `pip install -e .` |
| `cp2k-lsp-enhanced` | `poetry run pytest tests/` | `poetry install` |
| `abacus-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `abinit-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `gpumd-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `gromacs-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `lammps-lsp` | `cargo test` | `cargo build --release` |
| `pyscf-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `pyatb-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `mlip-lsp` | `python -m pytest tests/ -v` | `pip install -e .` |
| `cif-lsp` | `npm run test` | `npm run compile` |

## Supported File Types

| LSP Server | File Types |
|------------|-----------|
| `orca-lsp` | `.inp` (ORCA input) |
| `gamess-lsp` | `.inp` (GAMESS input) |
| `gaussian-lsp` | `.gjf`, `.com` (Gaussian input) |
| `nwchem-lsp` | `.nw` (NWChem input) |
| `qe-lsp` | `.in` (Quantum ESPRESSO input) |
| `VASP-LSP` | `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR` |
| `cp2k-lsp-enhanced` | `.inp` (CP2K input) |
| `abacus-lsp` | `INPUT`, `STRU`, `KPT` (ABACUS input) |
| `abinit-lsp` | `.abi`, `.abinit` (ABINIT input) |
| `gpumd-lsp` | `run.in`, `nep.in` (GPUMD input) |
| `gromacs-lsp` | `.top`, `.itp`, `.mdp`, `.gro` (GROMACS files) |
| `lammps-lsp` | `.lmp`, `.lammps`, `.lmps`, `in.lammps` (LAMMPS input) |
| `pyscf-lsp` | `.pyscf.py`, `run_pyscf.py` (PySCF scripts) |
| `pyatb-lsp` | `.pyatb.py`, `run_pyatb.py` (PyATB scripts) |
| `mlip-lsp` | `.mlip.json`, `.mlip.yaml`, `.mlip.yml` (MLIP config) |
| `cif-lsp` | `.cif` (Crystallographic Information File) |

## MatMaster Gap Servers

The following 9 LSP servers fill coverage gaps for materials-science and computational-physics workflows. They are registered in the bundled LSP registry (`src/lsp/registry.ts`) and are all marked **experimental** pending OpenQC integration testing.

### Server Details

| Server | Executable | Language ID | File Patterns | Install Source | Stability | Local Launch | Repository |
|--------|-----------|-------------|---------------|----------------|-----------|--------------|------------|
| `abacus-lsp` | `abacus-lsp` | `abacus` | `INPUT`, `STRU`, `KPT` | `pip install -e .` | Experimental | `pythonFunction` (`abacus_lsp.cli:lsp_main`) | [newtontech/abacus-lsp](https://github.com/newtontech/abacus-lsp) |
| `abinit-lsp` | `abinit-lsp` | `abinit` | `*.abi`, `*.abinit` | `pip install -e .` | Experimental | `pythonFunction` (`abinit_lsp.cli:lsp_main`) | [newtontech/abinit-lsp](https://github.com/newtontech/abinit-lsp) |
| `gpumd-lsp` | `gpumd-lsp` | `gpumd` | `run.in`, `nep.in` | `pip install -e .` | Experimental | `pythonFunction` (`gpumd_lsp.cli:lsp_main`) | [newtontech/gpumd-lsp](https://github.com/newtontech/gpumd-lsp) |
| `gromacs-lsp` | `gromacs-lsp` | `gromacs` | `*.top`, `*.itp`, `*.mdp`, `*.gro` | `pip install -e .` | Experimental | `pythonFunction` (`gromacs_lsp.cli:lsp_main`) | [newtontech/gromacs-lsp](https://github.com/newtontech/gromacs-lsp) |
| `lammps-lsp` | `lmp-lsp` | `lammps` | `*.lmp`, `*.lammps`, `*.lmps`, `in.lammps` | `cargo build --release` | Experimental | `cargoBinary` (`lmp-lsp`) | [newtontech/lammps-lsp](https://github.com/newtontech/lammps-lsp) |
| `pyscf-lsp` | `pyscf-lsp` | `pyscf` | `*.pyscf.py`, `run_pyscf.py` | `pip install -e .` | Experimental | `pythonFunction` (`pyscf_lsp.cli:lsp_main`) | [newtontech/pyscf-lsp](https://github.com/newtontech/pyscf-lsp) |
| `pyatb-lsp` | `pyatb-lsp` | `pyatb` | `*.pyatb.py`, `run_pyatb.py` | `pip install -e .` | Experimental | `pythonFunction` (`pyatb_lsp.cli:lsp_main`) | [newtontech/pyatb-lsp](https://github.com/newtontech/pyatb-lsp) |
| `mlip-lsp` | `mlip-lsp` | `mlip` | `*.mlip.json`, `*.mlip.yaml`, `*.mlip.yml` | `pip install -e .` | Experimental | `pythonFunction` (`mlip_lsp.cli:lsp_main`) | [newtontech/mlip-lsp](https://github.com/newtontech/mlip-lsp) |
| `cif-lsp` | `cif-lsp` | `cif` | `*.cif` | `npm install && npm run compile` | Experimental | `nodeScript` (`server/out/server.js`) | [newtontech/cif-lsp](https://github.com/newtontech/cif-lsp) |

### OpenQC Integration Issues

Each gap server has a corresponding tracking issue for full OpenQC integration:

| Server | Tracking Issue |
|--------|---------------|
| `abacus-lsp` | [newtontech/abacus-lsp#1](https://github.com/newtontech/abacus-lsp/issues/1) |
| `abinit-lsp` | [newtontech/abinit-lsp#1](https://github.com/newtontech/abinit-lsp/issues/1) |
| `gpumd-lsp` | [newtontech/gpumd-lsp#1](https://github.com/newtontech/gpumd-lsp/issues/1) |
| `gromacs-lsp` | [newtontech/gromacs-lsp#1](https://github.com/newtontech/gromacs-lsp/issues/1) |
| `lammps-lsp` | [newtontech/lammps-lsp#1](https://github.com/newtontech/lammps-lsp/issues/1) |
| `pyscf-lsp` | [newtontech/pyscf-lsp#1](https://github.com/newtontech/pyscf-lsp/issues/1) |
| `pyatb-lsp` | [newtontech/pyatb-lsp#1](https://github.com/newtontech/pyatb-lsp/issues/1) |
| `mlip-lsp` | [newtontech/mlip-lsp#1](https://github.com/newtontech/mlip-lsp/issues/1) |
| `cif-lsp` | [newtontech/cif-lsp#1](https://github.com/newtontech/cif-lsp/issues/1) |

## Shared Diagnostic JSON Contract

All LSP servers in the newtontech family produce diagnostics through a shared contract. OpenQC consumes these via the `DiagnosticsProvider` (`src/providers/lsp/DiagnosticsProvider.ts`), which converts them into VS Code `Diagnostic` objects.

### Parser Validation Result

Each parser (extending `BaseParser` in `src/parsers/base.ts`) returns a `ValidationResult`:

```typescript
interface ParseError {
  message: string;
  line: number;
  severity: 'error' | 'warning';
}

interface ParseWarning {
  message: string;
  line: number;
}

interface ValidationResult {
  valid: boolean;
  errors: ParseError[];
  warnings: ParseWarning[];
}
```

### Diagnostic Mapping

The `DiagnosticsProvider` maps each entry to a VS Code `Diagnostic`:

- `ParseError` entries with `severity: 'error'` become `DiagnosticSeverity.Error`
- `ParseWarning` entries become `DiagnosticSeverity.Warning`
- `source` is set to `"OpenQC"`
- `code` is set to `"<languageId>-error"` or `"<languageId>-warning"`

### Gap Server Conformance

All 9 gap servers follow the same `pythonFunction` launch contract (except `lammps-lsp` which uses `cargoBinary` and `cif-lsp` which uses `nodeScript`). Their diagnostic output must conform to the `ValidationResult` shape above to integrate correctly with the OpenQC `DiagnosticsProvider`.

## Alignment Policy

When changing diagnostics, completions, hover text, file detection, parser fixtures, or formatting in any LSP repo, also update or open an alignment issue in `OpenQC-VSCode` so the extension behavior stays consistent across all LSP servers.

---

Last updated: 2026-06-10
