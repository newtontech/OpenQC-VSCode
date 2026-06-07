# OpenQC

<div align="center">

**Your Universal Gateway to Computational Chemistry**

*Parse, visualize, and analyze ANY DFT/quantum chemistry/MD input files with one powerful extension*

[![Install](https://img.shields.io/badge/VS%20Code-Install-blue?logo=visual-studio-code)](https://marketplace.visualstudio.com/items?itemName=newtontech.openqc)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

[Features](#-supported-software) • [Quick Start](#-quick-start) • [Gallery](#-gallery) • [Roadmap](#-whats-coming)

</div>

---

## 🧪 Why OpenQC?

Tired of switching between different editors for VASP, Gaussian, ORCA, CP2K, and other quantum chemistry software? **OpenQC** brings everything together in one place.

**OpenQC automatically recognizes and parses your computational chemistry files** — whether you're doing DFT calculations, molecular dynamics, or wavefunction analysis. Just open your file and start working.

---

## 🎯 Supported Software

OpenQC is the VS Code-facing workspace for the newtontech computational chemistry LSP family. It should stay aligned with the standalone language servers in `newtontech/gaussian-lsp`, `newtontech/orca-lsp`, `newtontech/gamess-lsp`, `newtontech/qe-lsp`, and `newtontech/cp2k-lsp-enhanced`.

The repository root is the canonical VS Code extension package root. Development, tests, packaging, and Marketplace publishing should use the root `package.json`; stale nested extension package roots are not maintained.

We support **10+ major computational chemistry packages** with more on the way:

### LSP Alignment Matrix

| Format | Standalone LSP | OpenQC role |
|--------|----------------|-------------|
| Gaussian | `newtontech/gaussian-lsp` | Extension integration, syntax, visualization entry points |
| ORCA | `newtontech/orca-lsp` | Extension integration, syntax, visualization entry points |
| GAMESS (US) | `newtontech/gamess-lsp` | Extension integration, snippets, validation entry points |
| Quantum ESPRESSO | `newtontech/qe-lsp` | Extension integration and shared UX conventions |
| CP2K | `newtontech/cp2k-lsp-enhanced` | Extension integration with CP2K parser/linter tooling |
| VASP / NWChem / others | OpenQC modules or future LSPs | Keep APIs compatible with the same LSP contract |

### ✅ Fully Supported (Now)

| Software | File Types | Features |
|----------|-----------|----------|
| **VASP** | `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR` | Visualization + Syntax |
| **Gaussian** | `.com`, `.gjf` | Visualization + Syntax |
| **ORCA** | `.inp` | Visualization + Syntax |
| **CP2K** | `.inp` | Visualization + Syntax |
| **Quantum ESPRESSO** | `.in`, `.pw.in`, `.relax.in` | Visualization + Syntax |
| **GAMESS (US)** | `.inp` | Visualization + Syntax |
| **NWChem** | `.nw`, `.nwinp` | Visualization + Syntax |
| **Q-Chem** | `.in` | Syntax Highlighting |
| **ADF** | `.adf`, `.adfinput` | Syntax Highlighting |
| **TeraChem** | `.inp` | Syntax Highlighting |

### 🚧 Coming Soon

- **Molpro** - High-accuracy quantum chemistry
- **Psi4** - Open-source quantum chemistry
- **Molcas/OpenMolcas** - Multiconfigurational methods
- **DALTON** - Molecular properties
- **Turbomole** - Efficient DFT calculations
- **Crystal** - Periodic systems
- **Castep** - Materials modeling
- **Abinit** - Package for pseudopotential calculations

---

## ✨ What Can You Do?

### 🔬 Visualize Molecules in 3D

Open any input file and instantly see your molecular structure in beautiful 3D:

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

- **Syntax highlighting** for all supported formats
- **Error detection** as you type
- **Parameter validation** before you run
- **Smart completion** for common keywords

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
- [ ] Demand-proven ecosystem integrations from the [roadmap](PLAN.md)

## Marketplace Release Checklist

- Confirm `package.json` version, publisher, display name, icon, and keywords.
- Run the TypeScript build and extension packaging command from a clean checkout.
- Smoke test each aligned LSP family member: Gaussian, ORCA, GAMESS, QE, and CP2K.
- Smoke test VASP and NWChem sample files that are handled directly by OpenQC.
- Capture or refresh screenshots for syntax highlighting, diagnostics, 3D visualization, and validation.
- Publish release notes that list supported formats, known parser gaps, VS Code version tested, and LSP versions or commit SHAs.

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
