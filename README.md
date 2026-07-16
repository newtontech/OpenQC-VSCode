# OpenQC

<div align="center">

**VS Code integration for computational chemistry LSPs**

*Syntax, file detection, visualization entry points, and LSP startup for 17 computational chemistry and molecular-simulation formats*

Release identity: `newtontech.openqc@0.0.1` (`OpenQC - DFT/MD/Quantum Chemistry Suite`)

[![Install](https://img.shields.io/badge/VS%20Code-Install-blue?logo=visual-studio-code)](https://marketplace.visualstudio.com/items?itemName=newtontech.openqc)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

[Features](#-supported-software) • [Quick Start](#-quick-start) • [Gallery](#-gallery) • [Roadmap](#-whats-coming)

</div>

---

## 🧪 Why OpenQC?

Tired of switching between different editors for VASP, Gaussian, ORCA, CP2K, and other quantum chemistry software? **OpenQC** brings a shared VS Code surface to the newtontech language-server family.

**OpenQC recognizes supported computational chemistry and molecular-simulation files** and routes them to the matching language contribution, parser, viewer entry point, or configured LSP command.

---

## 🎯 Supported Software

OpenQC is the VS Code-facing workspace for the newtontech computational chemistry and molecular-simulation LSP family. It should stay aligned with the standalone language servers listed in `src/lsp/registry.ts` and summarized in `docs/LSP_COMPATIBILITY.md`.

The repository root is the canonical VS Code extension package root. Development, tests, packaging, and Marketplace publishing should use the root `package.json`; stale nested extension package roots are not maintained.

OpenQC currently wires **17 bundled LSP integrations** for computational chemistry, materials, and molecular-simulation workflows:

### Latest LSP Status

OpenQC tracks latest LSP support by the upstream default branch recorded in `src/lsp/registry.ts`. Run `npm run lsp:check-latest` before release or PR handoff to compare isolated source checkouts with configured remote branch heads. Run `npm run lsp:check-installed -- --strict` on a release host to additionally require the managed `.lsp-latest` runtime, server executable, agent CLI, source commit, manifest version, release-matrix version, and installed commit to agree. The installed-fleet check can emit the machine-readable `openqc.lsp.runtime-ledger.v1` report with `--json` or `--report-path <path>`.

### LSP Alignment Matrix

> For per-server parser status, diagnostics, completion, hover, formatting, code actions, and build commands, see the **[LSP Compatibility Matrix](docs/LSP_COMPATIBILITY.md)**.

| Format | Standalone LSP | OpenQC role |
|--------|----------------|-------------|
| ABACUS | `newtontech/abacus-lsp` | Language contribution, syntax, file detection, LSP startup |
| ABINIT | `newtontech/abinit-lsp` | Language contribution, syntax, file detection, LSP startup |
| CIF | `newtontech/cif-lsp` | Language contribution, syntax, file detection, LSP startup |
| CP2K | `newtontech/cp2k-lsp-enhanced` | Language contribution, syntax, file detection, LSP startup |
| VASP | `newtontech/VASP-LSP` | Language contribution, syntax, file detection, LSP startup |
| Gaussian | `newtontech/gaussian-lsp` | Language contribution, syntax, file detection, LSP startup |
| ORCA | `newtontech/orca-lsp` | Language contribution, syntax, file detection, LSP startup |
| GAMESS (US) | `newtontech/gamess-lsp` | Language contribution, syntax, file detection, LSP startup |
| Quantum ESPRESSO | `newtontech/qe-lsp` | Language contribution, syntax, file detection, LSP startup |
| NWChem | `newtontech/nwchem-lsp` | Language contribution, syntax, file detection, LSP startup |
| GPUMD | `newtontech/gpumd-lsp` | Language contribution, syntax, file detection, LSP startup |
| GROMACS | `newtontech/gromacs-lsp` | Language contribution, syntax, file detection, LSP startup |
| LAMMPS | `newtontech/lammps-lsp` | Language contribution, syntax, file detection, LSP startup |
| MLIP | `newtontech/mlip-lsp` | Language contribution, syntax, file detection, LSP startup |
| PyATB | `newtontech/pyatb-lsp` | Language contribution, syntax, file detection, LSP startup |
| PySCF | `newtontech/pyscf-lsp` | Language contribution, syntax, file detection, LSP startup |
| DP-GEN | `newtontech/dpgen-lsp` | Language contribution, syntax, file detection, LSP startup |

### Supported Integrations

| Software | File Types | Features |
|----------|-----------|----------|
| **ABACUS** | `INPUT`, `STRU`, `KPT` | LSP startup + syntax + file detection |
| **ABINIT** | `.abi`, `.abinit` | LSP startup + syntax + file detection |
| **CIF** | `.cif` | LSP startup + syntax + file detection |
| **CP2K** | `.inp` | LSP startup + syntax + file detection |
| **VASP** | `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR` | LSP startup + syntax + file detection |
| **Gaussian** | `.com`, `.gjf` | LSP startup + syntax + file detection |
| **ORCA** | `.inp` | LSP startup + syntax + file detection |
| **Quantum ESPRESSO** | `.in`, `.pw.in`, `.relax.in` | LSP startup + syntax + file detection |
| **GAMESS (US)** | `.inp` | LSP startup + syntax + file detection |
| **NWChem** | `.nw`, `.nwinp` | LSP startup + syntax + file detection |
| **GPUMD** | `run.in`, `nep.in` | LSP startup + syntax + file detection |
| **GROMACS** | `.top`, `.itp`, `.mdp`, `.gro` | LSP startup + syntax + file detection |
| **LAMMPS** | `.lmp`, `.lammps`, `.lmps` | LSP startup + syntax + file detection |
| **MLIP** | `.mlip.json`, `.mlip.yaml`, `.mlip.yml` | LSP startup + syntax + file detection |
| **PyATB** | `.pyatb.py`, `run_pyatb.py` | LSP startup + syntax + file detection |
| **PySCF** | `.pyscf.py`, `run_pyscf.py` | LSP startup + syntax + file detection |
| **DP-GEN** | `param.json`, `machine.json` | LSP startup + syntax + file detection |

### 🚧 Coming Soon

- **Molpro** - High-accuracy quantum chemistry
- **Psi4** - Open-source quantum chemistry
- **Molcas/OpenMolcas** - Multiconfigurational methods
- **DALTON** - Molecular properties
- **Turbomole** - Efficient DFT calculations
- **Crystal** - Periodic systems
- **Castep** - Materials modeling

---

## ✨ What Can You Do?

### 🔬 Visualize Molecules in 3D

Open supported molecule or structure-bearing files and use the viewer entry points to inspect coordinates in 3D:

- **Rotate, zoom, and pan** to explore your system
- **Multiple rendering styles** — ball-and-stick, space-filling, wireframe
- **Real-time preview** — see changes as you edit
- Support for **molecules, crystals, and surfaces**

### 📊 Analyze Your Calculations

Extract and visualize key data from your output files:

- SCF convergence plots
- Energy optimization progress
- Geometry optimization trajectories
- Molecular dynamics evolution

### 📝 Write Better Input Files

- **Syntax highlighting** for registered formats
- **Diagnostics** from the configured standalone LSP when that server provides them
- **Parameter validation** where the corresponding parser or LSP implements it
- **Completion and hover** where the corresponding standalone LSP implements them

### 🗂️ Organize Your Work

- Built-in **Molecules sidebar** to track your systems
- **Job tracking** panel for monitoring calculations
- Quick access to recent files and projects

---

## 🚀 Quick Start

### 1. Install OpenQC

Search for **"OpenQC"** in the VS Code Extensions panel and click Install.

### 2. Open Your File

Open any computational chemistry file:

```
POSCAR          # VASP structure
job.com         # Gaussian input
calc.inp        # ORCA/CP2K input
```

OpenQC **automatically detects** the file type and activates the right tools.

### 3. Visualize

Click the 🧪 icon in the editor toolbar to see your structure in 3D!

### 4. Analyze

Use the 📊 icon to plot your calculation data.

---

## 🎨 Gallery

### Molecular Visualization
*See your molecules come to life with interactive 3D rendering*

```
┌─────────────────────────────────────┐
│         [3D Viewer Panel]           │
│                                     │
│      ◯────◯────◯                   │
│     ╱      ╲     ╲                  │
│    ◯        ◯───◯                  │
│                                     │
│   Benzene • C₆H₆ • 12 atoms         │
└─────────────────────────────────────┘
```

### Syntax Highlighting
*Your input files, beautifully formatted*

```
&FORCE_EVAL
  SUBSYS
    &KIND O
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND
    &COORD
      O  0.000000  0.000000  0.000000
      H  0.758602  0.000000 -0.504284
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
```

---

## 🔧 Configuration

OpenQC works out of the box, but you can customize it:

```json
{
  // Auto-open visualization when opening files
  "openqc.visualization.autoOpen": true,

  // Your preferred rendering engine
  "openqc.visualization.moleculeRenderer": "3Dmol.js",

  // Auto-refresh sidebar views
  "openqc.sidebar.autoRefresh": true
}
```

---

## 💡 Use Cases

### For Computational Chemists
- **Prepare inputs** faster with syntax highlighting and validation
- **Visualize structures** before submitting jobs
- **Debug convergence** issues with interactive plots

### For Experimentalists
- **Inspect computational models** shared by collaborators
- **Understand output** from quantum chemistry calculations
- **Prepare structures** for computational studies

### For Students & Educators
- **Learn quantum chemistry** with visual feedback
- **Understand input formats** with syntax highlighting
- **Explore molecular systems** interactively

### For Software Developers
- **Build tools** on top of OpenQC's parsing capabilities
- **Integrate** with your computational workflows
- **Extend** support for additional software

---

## 🌟 What's Coming?

### Near Term (v2.1)
- [ ] Format conversion between different quantum chemistry formats
- [ ] Batch processing — visualize multiple structures at once
- [ ] Export high-resolution images for publications
- [ ] Custom color schemes and rendering options

### Medium Term (v2.5)
- [ ] Real-time calculation monitoring
- [ ] Remote workflow documentation using existing VS Code capabilities
- [ ] Parameter templates and wizards
- [ ] Community examples for common calculation workflows

### Long Term (v3.0)
- [ ] AI-powered parameter optimization
- [ ] Natural language input generation
- [ ] Workflow automation
- [ ] Demand-proven ecosystem integrations from the [roadmap](docs/project/PLAN.md)

## Marketplace Release Checklist

- Confirm `package.json` version, publisher, display name, icon, and keywords.
- Run the TypeScript build and extension packaging command from a clean checkout.
- Run `npm run lsp:check-latest` and keep every bundled LSP at the configured remote branch head or an isolated latest worktree.
- Run `npm run lsp:check-installed -- --strict` after rebuilding `.lsp-latest`; a latest source checkout does not prove the executable on `PATH` came from that commit.
- Smoke test each bundled LSP integration listed in `docs/LSP_COMPATIBILITY.md`.
- Capture or refresh screenshots for syntax highlighting, diagnostics, 3D visualization, and validation.
- Publish release notes that list supported formats, known parser gaps, VS Code version tested, and LSP versions or commit SHAs.
- Use Node 22 and the locked dependencies: `nvm use && npm ci`.
- Run `make format lint typecheck test check`; `npm run check:release` additionally creates and verifies `vsix/openqc-0.0.1.vsix`, its SHA-256 file, and CycloneDX SBOM.
- Create `v0.0.1` only from a clean commit exactly synchronized with `origin/master`; `npm run release:tag` enforces that boundary.
- A pushed `v*` tag runs build-and-verify first, pauses at the protected `marketplace` environment, publishes the already verified VSIX, and creates the GitHub Release only after Marketplace publication succeeds.
- Never publish with a global `vsce`; the repository pins `@vscode/vsce` and the workflow invokes that local binary.

See [Marketplace Publishing](docs/MARKETPLACE-PUBLISHING.md) for the operator runbook. Creating or pushing a tag is a release action and is not part of ordinary PR verification.

## Issue Triage Policy

Use labels to keep the roadmap readable:

- `lsp-alignment`: behavior that must stay consistent with standalone newtontech LSP repositories.
- `parser`: file parsing, validation, or format detection.
- `visualization`: 3D viewer, rendering, structures, or export.
- `language-support`: syntax highlighting, completion, diagnostics.
- `marketplace`: packaging, install, icon, metadata, release notes.
- `good-first-issue`: small parser fixtures, docs, examples, or screenshots.

Every bug should include a minimal input file or a redacted snippet that reproduces the behavior.

---

## 🤝 Contributing

We welcome contributions! See something missing? Let us know:

- **Add support for your favorite quantum chemistry software**
- **Improve parsing for existing formats**
- **Enhance visualization features**
- **Fix bugs and improve performance**

[Contributing Guidelines →](CONTRIBUTING.md)

---

## 📚 Resources

- [Documentation](https://docs.openqc.dev)
- [API Reference](https://api.openqc.dev)
- [Report Issues](https://github.com/newtontech/OpenQC-VSCode/issues)
- [Feature Requests](https://github.com/newtontech/OpenQC-VSCode/discussions)

---

## 📄 Citation

If OpenQC helps your research, please cite us:

```bibtex
@software{openqc2026,
  title = {OpenQC: Universal VS Code Extension for Computational Chemistry},
  author = {NewtonTech},
  year = {2026},
  version = {2.0},
  url = {https://github.com/newtontech/OpenQC-VSCode}
}
```

---

## 📜 License

MIT License — see [LICENSE](LICENSE) for details.

---

<div align="center">

**Made with ❤️ for the computational chemistry community**

[⭐ Star us on GitHub](https://github.com/newtontech/OpenQC-VSCode) •
[💬 Join the discussion](https://github.com/newtontech/OpenQC-VSCode/discussions)

</div>
