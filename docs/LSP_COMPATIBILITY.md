# LSP Compatibility Matrix

This document is generated from `src/lsp/registry.ts` and the LSP family gate report.
Regenerate it with `npm run lsp:generate-compatibility-doc` after registry, backend, or gate changes.

## Release Gate Summary

| Metric | Value |
|--------|-------|
| Backends | 17 |
| Passing | 17 |
| Blocking gaps | 0 |
| Warnings | 0 |
| Graduation score | 100/100 |

## Backend Matrix

| Backend | Language | Branch | File Types | Registry Stability | Gate Maturity | Release Evidence | Manifest | Provenance | Fixtures | Smoke | Diagnostics | Traceability |
|---------|----------|--------|------------|--------------------|---------------|------------------|----------|------------|----------|-------|-------------|--------------|
| [newtontech/abacus-lsp](https://github.com/newtontech/abacus-lsp) | `abacus` | `main` | `INPUT`, `STRU`, `KPT` | experimental | stable | `v0.1.0`<br>VERSION `0.1.0`<br>HEAD `71e0ead0463d` | pass | pass | pass | pass | pass | pass |
| [newtontech/abinit-lsp](https://github.com/newtontech/abinit-lsp) | `abinit` | `main` | `.abi`, `.abinit` | experimental | stable | `v0.1.0`<br>VERSION `0.1.0`<br>HEAD `b1c3124534ec` | pass | pass | pass | pass | pass | pass |
| [newtontech/cif-lsp](https://github.com/newtontech/cif-lsp) | `cif` | `master` | `.cif` | experimental | stable | `v1.0.2`<br>HEAD `314942468eff` | pass | pass | pass | pass | pass | pass |
| [newtontech/cp2k-lsp-enhanced](https://github.com/newtontech/cp2k-lsp-enhanced) | `cp2k` | `develop` | `.inp` | experimental | stable | `v0.9.1`<br>HEAD `7661a8c277fa` | pass | pass | pass | pass | pass | pass |
| [newtontech/VASP-LSP](https://github.com/newtontech/VASP-LSP) | `vasp` | `main` | `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR`, `CONTCAR`, `OSZICAR`, `OUTCAR`, `vasprun.xml` | stable | stable | `v0.4.4`<br>VERSION `0.4.4`<br>HEAD `be948e700924` | pass | pass | pass | pass | pass | pass |
| [newtontech/gaussian-lsp](https://github.com/newtontech/gaussian-lsp) | `gaussian` | `main` | `.gjf`, `.com` | stable | stable | `v0.2.11`<br>VERSION `0.2.11`<br>HEAD `282fb98ead89` | pass | pass | pass | pass | pass | pass |
| [newtontech/orca-lsp](https://github.com/newtontech/orca-lsp) | `orca` | `main` | `.inp` | stable | stable | `v0.5.5`<br>VERSION `0.5.5`<br>HEAD `7e7fa2f5d7ca` | pass | pass | pass | pass | pass | pass |
| [newtontech/qe-lsp](https://github.com/newtontech/qe-lsp) | `qe` | `main` | `.in`, `.pw.in`, `.relax.in`, `.vc-relax.in`, `.scf.in`, `.nscf.in`, `.bands.in`, `.ph.in`, `.dos.in` | stable | stable | `v0.1.0`<br>VERSION `0.1.0`<br>HEAD `c54045ae3b42` | pass | pass | pass | pass | pass | pass |
| [newtontech/gamess-lsp](https://github.com/newtontech/gamess-lsp) | `gamess` | `main` | `.inp` | stable | stable | `v0.1.0`<br>VERSION `0.1.0`<br>HEAD `289f1ac125fd` | pass | pass | pass | pass | pass | pass |
| [newtontech/nwchem-lsp](https://github.com/newtontech/nwchem-lsp) | `nwchem` | `main` | `.nw`, `.nwinp` | experimental | stable | `v0.3.0`<br>HEAD `cecc5ac62fb2` | pass | pass | pass | pass | pass | pass |
| [newtontech/gpumd-lsp](https://github.com/newtontech/gpumd-lsp) | `gpumd` | `main` | `run.in`, `nep.in` | experimental | stable | `v0.1.0`<br>VERSION `0.1.0`<br>HEAD `6c5d623c9d79` | pass | pass | pass | pass | pass | pass |
| [newtontech/gromacs-lsp](https://github.com/newtontech/gromacs-lsp) | `gromacs` | `main` | `.top`, `.itp`, `.mdp`, `.gro` | experimental | stable | `v0.0.3`<br>VERSION `0.0.3`<br>HEAD `61623cc06e71` | pass | pass | pass | pass | pass | pass |
| [newtontech/lammps-lsp](https://github.com/newtontech/lammps-lsp) | `lammps` | `master` | `.lmp`, `.lammps`, `.lmps`, `in.lammps` | experimental | stable | `0.1.0-pre-release-3`<br>HEAD `6f085a20b8f7` | pass | pass | pass | pass | pass | pass |
| [newtontech/mlip-lsp](https://github.com/newtontech/mlip-lsp) | `mlip` | `main` | `.mlip.json`, `.mlip.yaml`, `.mlip.yml`, `mlip.json`, `mlip.yaml`, `mlip.yml` | experimental | stable | `v0.2.0`<br>VERSION `0.2.0`<br>HEAD `d428ba5fbbda` | pass | pass | pass | pass | pass | pass |
| [newtontech/pyatb-lsp](https://github.com/newtontech/pyatb-lsp) | `pyatb` | `main` | `.pyatb.py`, `run_pyatb.py` | experimental | stable | `v0.1.0`<br>VERSION `0.1.0`<br>HEAD `b100e3ef8f12` | pass | pass | pass | pass | pass | pass |
| [newtontech/pyscf-lsp](https://github.com/newtontech/pyscf-lsp) | `pyscf` | `main` | `.pyscf.py`, `run_pyscf.py` | experimental | stable | `v0.1.0`<br>VERSION `0.1.0`<br>HEAD `fe377735407b` | pass | pass | pass | pass | pass | pass |
| [newtontech/dpgen-lsp](https://github.com/newtontech/dpgen-lsp) | `dpgen` | `main` | `param.json`, `machine.json` | experimental | stable | `v0.0.1`<br>VERSION `0.1.0`<br>HEAD `436419c1ced1` | pass | pass | pass | pass | pass | pass |

## Agent CLI Readiness

| Backend | Agent CLI | Operations | Help Smoke | Closed Loop |
|---------|-----------|------------|------------|-------------|
| `abacus-lsp` | `abacus-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `abinit-lsp` | `abinit-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `cif-lsp` | `cif-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `cp2k-lsp-enhanced` | `cp2k-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `vasp-lsp` | `vasp-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `gaussian-lsp` | `gaussian-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `orca-lsp` | `orca-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `qe-lsp` | `qe-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `gamess-lsp` | `gamess-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `nwchem-lsp` | `nwchem-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `gpumd-lsp` | `gpumd-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `gromacs-lsp` | `gromacs-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `lammps-lsp` | `lammps-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `mlip-lsp` | `mlip-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `pyatb-lsp` | `pyatb-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `pyscf-lsp` | `pyscf-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |
| `dpgen-lsp` | `dpgen-lsp-tool` | `check`, `context`, `complete`, `hover`, `symbols`, `fix` | pass | partial |

## Actionable Gate Gaps

No blocking or warning gaps are reported by `scripts/check-lsp-family.mjs`.

## OpenQC Integration Guarantees

| Capability | OpenQC guarantee |
|------------|------------------|
| Language contribution | Every registry entry has a matching `contributes.languages` item in `package.json`. |
| Grammar contribution | Every language has a TextMate grammar under `syntaxes/`. |
| Configuration | Every language exposes `enabled`, `path`, `command`, `args`, and `env` settings under `openqc.lsp.<languageId>.*`. |
| Startup | OpenQC starts the configured command over stdio and can prefer sibling local repositories when no user override is set. |
| Latest tracking | `npm run lsp:check-latest -- --fail-on-drift` verifies local checkout cleanliness, remote HEAD parity, and agent CLI help. |
| Family maturity | `npm run lsp:check-family -- --strict` verifies manifest, provenance, fixture, smoke, and Diagnostic Engine readiness. |
| Traceability aggregation | `npm run lsp:check-family -- --strict` consumes backend docstring/wiki/raw traceability reports when present and reports missing or broken evidence links. |
| Bohrium routing | `npm run lsp:check-bohrium-registry` verifies OpenQC and Bohrium backend ids and agent CLI names stay aligned. |

## Known Limits

- This gate proves OpenQC integration readiness and the shared agent-facing contract; it does not prove exhaustive grammar coverage for every historical software release.
- Complete version-aware parser/rule coverage and official-doc ingestion remain tracked in each backend repository issue queue.
- Runtime scientific correctness still requires backend-specific fixture suites and, where available, real executable/log smoke tests.
- OpenQC launches the configured executable or sibling checkout; it does not vendor or pin standalone backend binaries.
