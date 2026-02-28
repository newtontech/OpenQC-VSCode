# OpenQC-VSCode

> VSCode Extension for Quantum Chemistry Input Files: Visualization, Editing, Conversion, and AI-Assisted Modification

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Build Status](https://github.com/newtontech/OpenQC-VSCode/workflows/CI%2FCD%20Pipeline/badge.svg)](https://github.com/newtontech/OpenQC-VSCode/actions)
[![Coverage](https://codecov.io/gh/newtontech/OpenQC-VSCode/branch/main/graph/badge.svg)](https://codecov.io/gh/newtontech/OpenQC-VSCode)

## Features

### 🎯 Core Capabilities

- **Universal Parser**: Support for VASP, Gaussian, ORCA, and more
- **3D Visualization**: Interactive molecular structure rendering
- **Format Conversion**: Seamless conversion between different quantum chemistry formats
- **AI Assistance**: Smart suggestions and optimizations powered by AI
- **Syntax Highlighting**: Beautiful syntax highlighting for all supported formats

### 🔬 Supported Formats

| Format | Read | Write | Visualize | Convert |
|--------|------|-------|-----------|---------|
| VASP (INCAR, POSCAR, KPOINTS, POTCAR) | ✅ | ✅ | ✅ | ✅ |
| Gaussian (.com, .gjf) | ✅ | ✅ | ✅ | ✅ |
| ORCA (.inp) | ✅ | ✅ | ✅ | ✅ |
| XYZ | ✅ | ✅ | ✅ | ✅ |
| PDB | ✅ | ✅ | ✅ | ✅ |
| CIF | ✅ | ✅ | ✅ | ✅ |

### 🤖 AI Features

- **Parameter Suggestions**: Get AI-powered recommendations for calculation parameters
- **Error Debugging**: Understand and fix failed calculations
- **Natural Language Generation**: Create input files from descriptions
- **Parameter Explanation**: Understand what each parameter does

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

### 2. Parse and Validate

The extension automatically parses your file and shows:
- Syntax highlighting
- Error markers for invalid syntax
- Hover tooltips for parameters

### 3. Visualize Structure

Press `Ctrl+Shift+P` and run `OpenQC: Visualize Structure` to see 3D rendering.

### 4. Convert Format

Run `OpenQC: Convert File` to convert between formats.

### 5. Get AI Suggestions

Run `OpenQC: AI Suggest` to get intelligent parameter recommendations.

## Configuration

Add these settings to your `settings.json`:

```json
{
  "openqc.ai.provider": "openai",
  "openqc.ai.model": "gpt-4",
  "openqc.ai.apiKey": "your-api-key",
  "openqc.visualization.renderer": "ngl",
  "openqc.parsers.preferPython": false,
  "openqc.conversion.preserveMetadata": true
}
```

## Commands

| Command | Description |
|---------|-------------|
| `OpenQC: Parse File` | Parse current file |
| `OpenQC: Visualize Structure` | Open 3D visualization |
| `OpenQC: Convert File` | Convert to different format |
| `OpenQC: Validate File` | Check for errors |
| `OpenQC: AI Suggest` | Get AI suggestions |
| `OpenQC: AI Explain` | Explain parameters |

## Documentation

- [Project Roadmap](PLAN.md) - Development plan and milestones
- [TDD Guidelines](TDD-GUIDELINES.md) - Testing best practices
- [Task Management](TASK-MANAGEMENT.md) - How we organize work
- [Architecture](docs/architecture/ARCHITECTURE.md) - System design
- [API Reference](docs/api/API.md) - Developer documentation

## Development

### Prerequisites

- Node.js 18+
- Python 3.8+ (for backend integrations)
- VSCode 1.85+

### Setup

```bash
# Clone repository
git clone https://github.com/newtontech/OpenQC-VSCode.git
cd OpenQC-VSCode

# Install dependencies
npm install
pip install -e .[dev]

# Run tests
npm test
pytest

# Build
npm run compile
```

### Testing

We follow Test-Driven Development (TDD):

```bash
# Run unit tests
npm run test:unit

# Run integration tests
npm run test:integration

# Run e2e tests
npm run test:e2e

# Generate coverage
npm run test:coverage
```

See [TDD-GUIDELINES.md](TDD-GUIDELINES.md) for details.

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Ways to Contribute

- 🐛 Report bugs
- 💡 Suggest features
- 📝 Improve documentation
- 🔧 Submit pull requests
- ⭐ Star the repository

## Roadmap

See [PLAN.md](PLAN.md) for detailed roadmap.

### Current Phase: Foundation (Weeks 1-4)
- Core parsing engine
- Syntax highlighting
- Basic visualization

### Upcoming Phases
- Phase 2: Advanced visualization
- Phase 3: Format conversion
- Phase 4: AI integration
- Phase 5: Performance optimization

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- [dpdata](https://github.com/deepmodeling/dpdata) - Format conversion
- [ASE](https://wiki.fysik.dtu.dk/ase/) - Atomic Simulation Environment
- [NGL Viewer](http://nglviewer.org/) - 3D molecular visualization
- [Three.js](https://threejs.org/) - 3D graphics library

## Support

- 📧 Email: support@newtontech.com
- 💬 Discord: [Join our community](https://discord.gg/openqc)
- 📖 Documentation: [docs.openqc.dev](https://docs.openqc.dev)
- 🐛 Issues: [GitHub Issues](https://github.com/newtontech/OpenQC-VSCode/issues)

## Citation

If you use OpenQC-VSCode in your research, please cite:

```bibtex
@software{openqc2026,
  title = {OpenQC-VSCode: VSCode Extension for Quantum Chemistry Input Files},
  author = {NewtonTech},
  year = {2026},
  url = {https://github.com/newtontech/OpenQC-VSCode}
}
```

---

**Made with ❤️ by the NewtonTech team**
