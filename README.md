# OpenQC-VSCode

> **Universal VS Code extension for quantum chemistry software with multi-LSP support and advanced visualization**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Build Status](https://github.com/newtontech/OpenQC-VSCode/workflows/CI%2FCD%20Pipeline/badge.svg)](https://github.com/newtontech/OpenQC-VSCode/actions)
[![Coverage](https://codecov.io/gh/newtontech/OpenQC-VSCode/branch/main/graph/badge.svg)](https://codecov.io/gh/newtontech/OpenQC-VSCode/branch/main/graph/badge.svg)

## 🚀 What's New in v2.0

**OpenQC-VSCode is now a universal platform for quantum chemistry!** We've expanded from a single-tool focus to supporting **7 major quantum chemistry packages** with automatic LSP detection and unified visualization.

### Major Features
- ✅ **Universal LSP Support** - Auto-detect and launch language servers for 7 quantum chemistry packages
- ✅ **Molecular Visualization** - Interactive 3D rendering with 3Dmol.js
- ✅ **Data Visualization** - Plot SCF energies, convergence data with Plotly.js
- ✅ **Input Preview** - Structured preview of input file parameters
- ✅ **Developer Tools** - Syntax highlighting, validation, and LSP diagnostics

## Features

### 🎯 Core Capabilities

- **🚀 Universal LSP Support**: Automatic detection and management of language servers for:
  - CP2K
  - VASP
  - Gaussian
  - ORCA
  - Quantum ESPRESSO
  - GAMESS
  - NWChem

- **🔬 Molecular Visualization**: Interactive 3D structure rendering with 3Dmol.js
  - Multiple visualization styles: stick, sphere, line, cartoon
  - Spin and zoom controls
  - Real-time structure preview from input files

- **📊 Data Visualization**: Plot calculation data with Plotly.js
  - SCF energy convergence
  - K-point grids
  - Automatic data extraction from output files
  - Interactive and responsive charts

- **📝 Input Preview**: Structured display of input file parameters
  - Section-based organization
  - Parameter extraction and display
  - Syntax highlighting for all formats

- **🛠️ Developer Tools**:
  - Syntax highlighting for all 7 quantum chemistry formats
  - File type auto-detection
  - Language server management (start/stop/restart)
  - Error diagnostics from LSPs

### 🔬 Supported Quantum Chemistry Packages

| Package | Language ID | Files | Syntax | Visualization | LSP |
|---------|--------------|-------|--------|---------------|-----|
| **CP2K** | `cp2k` | `.inp` | ✅ | ✅ | ✅ |
| **VASP** | `vasp` | `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR` | ✅ | ✅ | ✅ |
| **Gaussian** | `gaussian` | `.gjf`, `.com` | ✅ | ✅ | ✅ |
| **ORCA** | `orca` | `.inp` | ✅ | ✅ | ✅ |
| **Quantum ESPRESSO** | `qe` | `.in`, `.pw.in`, `.relax.in`, etc. | ✅ | ✅ | ✅ |
| **GAMESS** | `gamess` | `.inp` | ✅ | ✅ | ✅ |
| **NWChem** | `nwchem` | `.nw`, `.nwinp` | ✅ | ✅ | ✅ |

## Installation

### From VSCode Marketplace

1. Open VSCode
2. Press `Ctrl+Shift+X` (Windows/Linux) or `Cmd+Shift+X` (macOS)
3. Search for "OpenQC-VSCode"
4. Click "Install"

### From Source

```bash
git clone https://github.com/newtontech/OpenQC-VSCode.git
cd OpenQC-VSCode
npm install
npm run compile
# Package and install locally
npx vsce package
code --install-extension openqc-vscode-*.vsix
```

## Quick Start

### 1. Open a Quantum Chemistry File

Open any supported file (e.g., `POSCAR`, `input.com`, `job.inp`) in VSCode.
The extension will **automatically detect** the file type and launch the appropriate language server.

### 2. Visualize Molecular Structure

- **Option 1**: Click the structure icon in the editor title bar
- **Option 2**: Press `Ctrl+Shift+P` and run `OpenQC: Visualize Structure`
- **Option 3**: Right-click and select "OpenQC: Visualize Structure"

### 3. Plot Calculation Data

- **Option 1**: Click the plot icon in the editor title bar
- **Option 2**: Press `Ctrl+Shift+P` and run `OpenQC: Plot Calculation Data`

### 4. Preview Input File

- **Option 1**: Press `Ctrl+Shift+P` and run `OpenQC: Preview Input File`
- **Option 2**: Right-click and select "OpenQC: Preview Input File"

### 5. Manage Language Servers

```bash
# Start language server
Ctrl+Shift+P > OpenQC: Start Language Server

# Stop language server
Ctrl+Shift+P > OpenQC: Stop Language Server

# Restart language server
Ctrl+Shift+P > OpenQC: Restart Language Server
```

## Configuration

### LSP Configuration

Configure language server paths in your `settings.json`:

```json
{
  "openqc.lsp.cp2k.enabled": true,
  "openqc.lsp.cp2k.path": "cp2k-lsp-enhanced",
  
  "openqc.lsp.vasp.enabled": true,
  "openqc.lsp.vasp.path": "vasp-lsp",
  
  "openqc.lsp.gaussian.enabled": true,
  "openqc.lsp.gaussian.path": "gaussian-lsp",
  
  "openqc.lsp.orca.enabled": true,
  "openqc.lsp.orca.path": "orca-lsp",
  
  "openqc.lsp.qe.enabled": true,
  "openqc.lsp.qe.path": "qe-lsp",
  
  "openqc.lsp.gamess.enabled": true,
  "openqc.lsp.gamess.path": "gamess-lsp",
  
  "openqc.lsp.nwchem.enabled": true,
  "openqc.lsp.nwchem.path": "nwchem-lsp"
}
```

### Visualization Configuration

```json
{
  "openqc.visualization.moleculeRenderer": "3Dmol.js",
  "openqc.visualization.plotLibrary": "Plotly.js",
  "openqc.visualization.autoOpen": true
}
```

## Commands

| Command | Description |
|---------|-------------|
| `OpenQC: Visualize Structure` | Open 3D molecular structure viewer |
| `OpenQC: Plot Calculation Data` | Plot SCF energies and convergence data |
| `OpenQC: Preview Input File` | Show structured preview of input file |
| `OpenQC: Start Language Server` | Manually start the language server |
| `OpenQC: Stop Language Server` | Stop the language server |
| `OpenQC: Restart Language Server` | Restart the language server |

## Architecture

### LSP Manager

The LSP Manager automatically:
1. Detects the quantum chemistry software from file extension and content
2. Launches the appropriate language server
3. Manages server lifecycle (start/stop/restart)
4. Handles multiple file types simultaneously

### Visualization Pipeline

```
Input File → Parser → Atoms/Data → Webview → 3Dmol.js/Plotly.js
```

### File Type Detection

Multi-layer detection:
1. **Filename match** - Exact filename (e.g., `INCAR`, `POSCAR`)
2. **Extension match** - File extension (e.g., `.inp`, `.gjf`)
3. **Content analysis** - Regex patterns for ambiguous cases

## Development

### Prerequisites

- Node.js 18+
- TypeScript 5.3+
- VSCode 1.85+

### Setup

```bash
# Clone repository
git clone https://github.com/newtontech/OpenQC-VSCode.git
cd OpenQC-VSCode

# Install dependencies
npm install

# Compile
npm run compile

# Watch mode for development
npm run watch
```

### Project Structure

```
OpenQC-VSCode/
├── src/
│   ├── extension.ts              # Extension entry point
│   ├── managers/
│   │   ├── LSPManager.ts         # Language server management
│   │   └── FileTypeDetector.ts   # File type detection
│   ├── providers/
│   │   ├── StructureViewer.ts    # 3D structure visualization
│   │   └── DataPlotter.ts        # Data plotting
│   └── visualizers/
│       └── Molecule3D.ts         # Molecule parsing
├── syntaxes/                      # Syntax highlighting
├── language-configurations/       # Language config
├── package.json
├── tsconfig.json
└── README.md
```

### Testing

```bash
# Run tests
npm test

# Run with coverage
npm run test:coverage

# Run linting
npm run lint
```

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Areas for Contribution

- 🐛 Bug fixes
- 💡 New features
- 📝 Documentation improvements
- 🎨 UI/UX improvements
- 🔧 Additional LSP integrations
- 📊 New visualization types

## Roadmap

### v2.1 (Planned)
- [ ] Format conversion between quantum chemistry formats
- [ ] Batch visualization
- [ ] Custom color schemes
- [ ] Export images

### v2.2 (Planned)
- [ ] Real-time calculation monitoring
- [ ] Integration with job schedulers
- [ ] Parameter templates
- [ ] Cloud storage integration

### v3.0 (Future)
- [ ] AI-powered parameter optimization
- [ ] Natural language input generation
- [ ] Workflow automation
- [ ] Multi-package job orchestration

## Documentation

- [Project Roadmap](PLAN.md) - Development plan and milestones
- [TDD Guidelines](TDD-GUIDELINES.md) - Testing best practices
- [Task Management](TASK-MANAGEMENT.md) - How we organize work
- [Architecture](docs/architecture/ARCHITECTURE.md) - System design
- [API Reference](docs/api/API.md) - Developer documentation

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **3Dmol.js** - Interactive molecular visualization
- **Plotly.js** - Interactive data visualization
- **vscode-languageclient** - VS Code Language Client library
- Quantum Chemistry Community - Inspiration and feedback

## Support

- 📧 Email: support@newtontech.com
- 💬 Discord: [Join our community](https://discord.gg/openqc)
- 📖 Documentation: [docs.openqc.dev](https://docs.openqc.dev)
- 🐛 Issues: [GitHub Issues](https://github.com/newtontech/OpenQC-VSCode/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/newtontech/OpenQC-VSCode/discussions)

## Citation

If you use OpenQC-VSCode in your research, please cite:

```bibtex
@software{openqc2026,
  title = {OpenQC-VSCode: Universal VS Code Extension for Quantum Chemistry},
  author = {NewtonTech},
  year = {2026},
  version = {2.0},
  url = {https://github.com/newtontech/OpenQC-VSCode}
}
```

---

**Made with ❤️ by the NewtonTech team**

⭐ **Star us on GitHub** if you find OpenQC-VSCode helpful!
