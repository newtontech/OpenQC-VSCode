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

## Alignment Policy

When changing diagnostics, completions, hover text, file detection, parser fixtures, or formatting in any LSP repo, also update or open an alignment issue in `OpenQC-VSCode` so the extension behavior stays consistent across all LSP servers.

---

Last updated: 2026-06-08
