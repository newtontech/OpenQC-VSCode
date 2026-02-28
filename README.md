# OpenQC

<div align="center">

**Open-Source VSCode Extension for Quantum Chemistry**

[![VSCode Version](https://img.shields.io/badge/VSCode-1.85+-blue.svg)](https://code.visualstudio.com/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub Stars](https://img.shields.io/github/stars/newtontech/OpenQC.svg)](https://github.com/newtontech/OpenQC/stargazers)

[Features](#features) • [Installation](#installation) • [Documentation](#documentation) • [Contributing](#contributing)

</div>

---

## Overview

OpenQC is a comprehensive VSCode extension designed to streamline quantum chemistry workflows. It provides visualization, editing, format conversion, and AI-assisted modification for all major quantum chemistry input file formats.

### Supported Formats

| Software | File Types | Status |
|----------|------------|--------|
| **Gaussian** | `.gjf`, `.com` | ✅ Full Support |
| **VASP** | `POSCAR`, `CONTCAR`, `INCAR` | ✅ Full Support |
| **Quantum ESPRESSO** | `.pw`, `.scf`, `.nscf` | ✅ Full Support |
| **ORCA** | `.inp` | ✅ Full Support |
| **CP2K** | `.inp` | ✅ Full Support |
| **Q-Chem** | `.qcinp` | 🔄 In Progress |
| **NWChem** | `.nw` | 🔄 In Progress |
| **Psi4** | `.dat` | 📋 Planned |

---

## Features

### 🔬 Visualization
- **3D Molecular Viewer**: Interactive visualization of molecular structures
- **Crystal Structure Visualization**: Unit cell and supercell rendering
- **Orbital Visualization**: Render molecular orbitals and electron density
- **Trajectory Animation**: Play MD trajectories directly in VSCode

### ✏️ Smart Editing
- **Syntax Highlighting**: Comprehensive syntax support for all QC formats
- **Auto-completion**: Intelligent suggestions for keywords and parameters
- **Parameter Validation**: Real-time checking for invalid parameters
- **Templates Library**: Quick-start templates for common calculations

### 🔄 Format Conversion
- **dpdata Integration**: Seamless conversion between 50+ formats
- **ASE Integration**: Atomic Simulation Environment support
- **Batch Conversion**: Convert multiple files at once
- **Structure Optimization**: Format-aware structure optimization

### 🤖 AI Assistance
- **Claude Code Integration**: Native MCP server for Claude Code
- **Copilot Enhancement**: AI-assisted parameter suggestions
- **Natural Language Queries**: "Change functional to B3LYP" works!
- **Error Analysis**: AI-powered error message interpretation

### 🖥️ Remote Computing
- **One-Click SSH**: Connect to compute servers instantly
- **Slurm Integration**: Submit and monitor jobs
- **Job Queue Management**: Track all your running jobs
- **File Transfer**: Drag-and-drop remote file access

### 🔌 AI Protocol Support
- **MCP Server**: Full Model Context Protocol implementation
- **ACP Protocol**: Agent Control Protocol for automation
- **REST API**: Programmatic access to all features
- **CLI Tools**: Command-line interface for scripting

---

## Installation

### From VSCode Marketplace
\`\`\`bash
# Coming soon to marketplace
code --install-extension newtontech.openqc
\`\`\`

### From Source
\`\`\`bash
# Clone the repository
git clone https://github.com/newtontech/OpenQC.git
cd OpenQC

# Install dependencies
cd vscode-extension
npm install

# Build and package
npm run compile
vsce package

# Install in VSCode
code --install-extension openqc-*.vsix
\`\`\`

### Quick Start
1. Open a quantum chemistry input file
2. Press \`Ctrl+Shift+P\` → "OpenQC: Visualize"
3. Explore your molecular structure!

---

## Project Structure

\`\`\`
OpenQC/
├── vscode-extension/     # VSCode extension (TypeScript)
│   ├── src/
│   │   ├── extension.ts      # Main entry point
│   │   ├── providers/        # Language providers
│   │   ├── views/            # Webview panels
│   │   └── commands/         # Command implementations
│   ├── syntaxes/             # TextMate grammars
│   ├── language-configuration/
│   └── package.json
│
├── core/                 # Core Python library
│   ├── openqc/
│   │   ├── parsers/          # File format parsers
│   │   ├── converters/       # Format converters
│   │   ├── visualizers/      # Visualization engine
│   │   └── utils/            # Utility functions
│   └── pyproject.toml
│
├── server/               # Remote compute server
│   ├── ssh_handler.py        # SSH connection management
│   ├── slurm_interface.py    # Slurm job management
│   └── file_transfer.py      # Remote file operations
│
├── ai-protocols/         # AI protocol implementations
│   ├── mcp-server/           # MCP server for Claude
│   │   ├── server.py
│   │   └── tools/
│   └── acp-adapter/          # ACP protocol adapter
│       └── adapter.py
│
├── docs/                 # Documentation
│   ├── user-guide/
│   ├── api-reference/
│   └── examples/
│
├── tests/                # Test suite
│   ├── unit/
│   ├── integration/
│   └── fixtures/
│
└── examples/             # Example input files
    ├── gaussian/
    ├── vasp/
    ├── qe/
    └── orca/
\`\`\`

---

## AI Integration

### Claude Code (MCP)

OpenQC provides a native MCP server for seamless Claude Code integration:

\`\`\`json
// Add to Claude Code config
{
  "mcpServers": {
    "openqc": {
      "command": "python",
      "args": ["/path/to/OpenQC/ai-protocols/mcp-server/server.py"]
    }
  }
}
\`\`\`

**Available MCP Tools:**
- \`parse_qc_file\` - Parse any QC input file
- \`convert_format\` - Convert between formats
- \`visualize_structure\` - Generate 3D visualization
- \`modify_parameters\` - AI-assisted parameter changes
- \`validate_input\` - Check for errors
- \`submit_job\` - Submit to remote servers

### Example Claude Code Session
\`\`\`
User: Parse this Gaussian file and convert it to VASP format
Claude: [Uses openqc.parse_qc_file, then openqc.convert_format]
Here's your VASP POSCAR file. The structure contains...

User: Change the functional to PBE0 and increase the basis set
Claude: [Uses openqc.modify_parameters]
Done! I've updated the functional to PBE0 and changed...
\`\`\`

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

<div align="center">

Made with ❤️ by the OpenQC Team

</div>
