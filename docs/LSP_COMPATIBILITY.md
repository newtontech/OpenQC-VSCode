# LSP Compatibility Matrix

This matrix is the OpenQC-facing contract for the newtontech LSP family. The source of truth is `src/lsp/registry.ts`; `package.json`, language configurations, TextMate grammars, and this document must stay aligned with that registry.

Last updated: 2026-06-15. Latest support is tracked against each repository default branch unless a release tag is listed below. OpenQC does not pin old server binaries; it launches the configured executable or sibling local repository, so users can run the newest server available in their environment.

## Quick Reference

| LSP Server                                                                      | Language ID | Default Branch | Latest Release/Local Version | Domain           | File Types                                                                                            | Stability    | OpenQC Integration                         |
| ------------------------------------------------------------------------------- | ----------- | -------------- | ---------------------------- | ---------------- | ----------------------------------------------------------------------------------------------------- | ------------ | ------------------------------------------ |
| [newtontech/abacus-lsp](https://github.com/newtontech/abacus-lsp)               | `abacus`    | `main`         | 0.1.0                        | ABACUS           | `INPUT`, `STRU`, `KPT`                                                                                | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/abinit-lsp](https://github.com/newtontech/abinit-lsp)               | `abinit`    | `main`         | 0.1.0                        | ABINIT           | `.abi`, `.abinit`                                                                                     | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/cif-lsp](https://github.com/newtontech/cif-lsp)                     | `cif`       | `master`       | v1.0.2                       | CIF              | `.cif`                                                                                                | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/cp2k-lsp-enhanced](https://github.com/newtontech/cp2k-lsp-enhanced) | `cp2k`      | `develop`      | v0.9.1                       | CP2K             | `.inp`                                                                                                | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/VASP-LSP](https://github.com/newtontech/VASP-LSP)                   | `vasp`      | `main`         | 0.4.4                        | VASP             | `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR`, `CONTCAR`, `OSZICAR`, `OUTCAR`, `vasprun.xml`                 | Stable       | Syntax, file detection, direct LSP startup |
| [newtontech/gaussian-lsp](https://github.com/newtontech/gaussian-lsp)           | `gaussian`  | `main`         | 0.2.11                       | Gaussian         | `.gjf`, `.com`                                                                                        | Stable       | Syntax, file detection, direct LSP startup |
| [newtontech/orca-lsp](https://github.com/newtontech/orca-lsp)                   | `orca`      | `main`         | 0.5.4                        | ORCA             | `.inp`                                                                                                | Stable       | Syntax, file detection, direct LSP startup |
| [newtontech/qe-lsp](https://github.com/newtontech/qe-lsp)                       | `qe`        | `main`         | 0.1.0                        | Quantum ESPRESSO | `.in`, `.pw.in`, `.relax.in`, `.vc-relax.in`, `.scf.in`, `.nscf.in`, `.bands.in`, `.ph.in`, `.dos.in` | Stable       | Syntax, file detection, direct LSP startup |
| [newtontech/gamess-lsp](https://github.com/newtontech/gamess-lsp)               | `gamess`    | `main`         | 0.1.0                        | GAMESS           | `.inp`                                                                                                | Stable       | Syntax, file detection, direct LSP startup |
| [newtontech/nwchem-lsp](https://github.com/newtontech/nwchem-lsp)               | `nwchem`    | `main`         | v0.3.0 / local 0.5.0         | NWChem           | `.nw`, `.nwinp`                                                                                       | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/gpumd-lsp](https://github.com/newtontech/gpumd-lsp)                 | `gpumd`     | `main`         | 0.1.0                        | GPUMD            | `run.in`, `nep.in`                                                                                    | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/gromacs-lsp](https://github.com/newtontech/gromacs-lsp)             | `gromacs`   | `main`         | v0.0.2                       | GROMACS          | `.top`, `.itp`, `.mdp`, `.gro`                                                                        | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/lammps-lsp](https://github.com/newtontech/lammps-lsp)               | `lammps`    | `master`       | 0.1.0-pre-release-3          | LAMMPS           | `.lmp`, `.lammps`, `.lmps`, `in.lammps`                                                               | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/mlip-lsp](https://github.com/newtontech/mlip-lsp)                   | `mlip`      | `main`         | 0.1.0                        | MLIP             | `.mlip.json`, `.mlip.yaml`, `.mlip.yml`, `mlip.json`, `mlip.yaml`, `mlip.yml`                         | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/pyatb-lsp](https://github.com/newtontech/pyatb-lsp)                 | `pyatb`     | `main`         | 0.1.0                        | PyATB            | `.pyatb.py`, `run_pyatb.py`                                                                           | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/pyscf-lsp](https://github.com/newtontech/pyscf-lsp)                 | `pyscf`     | `main`         | 0.1.0                        | PySCF            | `.pyscf.py`, `run_pyscf.py`                                                                           | Experimental | Syntax, file detection, direct LSP startup |
| [newtontech/dpgen-lsp](https://github.com/newtontech/dpgen-lsp)                 | `dpgen`     | `main`         | 0.1.0                        | DP-GEN           | `param.json`, `machine.json`                                                                          | Experimental | Syntax, file detection, direct LSP startup |

## Diagnostic Engine v1 Readiness

| LSP Server          | Diagnostic Engine | Agent CLI           | Rich Diagnostics        | Closed Loop |
| ------------------- | ----------------- | ------------------- | ----------------------- | ----------- |
| `abacus-lsp`        | v1                | `abacus-lsp-tool`   | Rich JSON check payload | Planned     |
| `abinit-lsp`        | v1                | `abinit-lsp-tool`   | Rich JSON check payload | Planned     |
| `cif-lsp`           | v1                | `cif-lsp-tool`      | Rich JSON check payload | Planned     |
| `cp2k-lsp-enhanced` | v1                | `cp2k-lsp-tool`     | Rich JSON check payload | Partial     |
| `dpgen-lsp`         | v1                | `dpgen-lsp-tool`    | Rich JSON check payload | Planned     |
| `gamess-lsp`        | v1                | `gamess-lsp-tool`   | Rich JSON check payload | Planned     |
| `gaussian-lsp`      | v1                | `gaussian-lsp-tool` | Rich JSON check payload | Planned     |
| `gpumd-lsp`         | v1                | `gpumd-lsp-tool`    | Rich JSON check payload | Planned     |
| `gromacs-lsp`       | v1                | `gromacs-lsp-tool`  | Rich JSON check payload | Planned     |
| `lammps-lsp`        | v1                | `lammps-lsp-tool`   | Rich JSON check payload | Partial     |
| `mlip-lsp`          | v1                | `mlip-lsp-tool`     | Rich JSON check payload | Planned     |
| `nwchem-lsp`        | v1                | `nwchem-lsp-tool`   | Rich JSON check payload | Planned     |
| `orca-lsp`          | v1                | `orca-lsp-tool`     | Rich JSON check payload | Planned     |
| `pyatb-lsp`         | v1                | `pyatb-lsp-tool`    | Rich JSON check payload | Planned     |
| `pyscf-lsp`         | v1                | `pyscf-lsp-tool`    | Rich JSON check payload | Planned     |
| `qe-lsp`            | v1                | `qe-lsp-tool`       | Rich JSON check payload | Partial     |
| `vasp-lsp`          | v1                | `vasp-lsp-tool`     | Rich JSON check payload | Partial     |

## Capability Manifest Contract

Every sibling LSP repository must expose `lsp-capabilities.json` at its repository root. OpenQC treats `src/lsp/registry.ts` as the bundled product default and treats each capability manifest as the current repository fact surface for agent-facing features.

The manifest records the registry id, language id, executable, default branch, file patterns, blocking policy, agent CLI operations, DiagnosticEnvelope schema path, LLM Wiki paths, smoke fixture paths, output-log patterns, and OpenQC context contract. The required agent CLI operations are `capabilities`, `check`, `context`, `complete`, `hover`, `symbols`, and `fix`, all with JSON output.

Run these checks before release:

```bash
npm run lsp:check-family
npm run lsp:check-latest
```

## OpenQC Integration Guarantees

| Capability            | OpenQC guarantee                                                                                                        |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Language contribution | Every registry entry has a matching `contributes.languages` item in `package.json`.                                     |
| Grammar contribution  | Every language has a TextMate grammar under `syntaxes/`.                                                                |
| Configuration         | Every language exposes `enabled`, `path`, `command`, `args`, and `env` settings under `openqc.lsp.<languageId>.*`.      |
| Startup               | OpenQC starts the configured command over stdio and can prefer sibling local repositories when no user override is set. |
| Latest tracking       | Every registry entry records the upstream default branch used for latest-version checks.                                |

## Release Verification

Before publishing a VSIX or Marketplace release:

1. Run `npm run typecheck`.
2. Run `npm run test:unit -- --runTestsByPath tests/unit/lsp/registry.test.ts`.
3. Run `npm run lsp:check-family` to verify every sibling manifest, agent CLI contract, wiki path, and fixture path.
4. Run `npm run lsp:check-latest` to compare sibling LSP checkouts with each configured remote default-branch HEAD.
5. Run `OpenQC LSP: Generate Compatibility Report` from VS Code and save the report with the release notes.
6. For each LSP, smoke test one valid file and one intentionally invalid file on the latest configured executable or sibling repository checkout.
7. Record the exact LSP version, release tag, or commit SHA used for the release smoke test.

## Known Limits

- The matrix records OpenQC integration support, not a promise that every standalone LSP implements every optional LSP method.
- `cp2k-lsp-enhanced` uses `develop` as its latest upstream branch.
- `cif-lsp` and `lammps-lsp` use `master` as their latest upstream branch.
- Experimental entries are wired for startup and language support, but may have narrower standalone server capabilities than the stable chemistry LSPs.
