# OpenQC VSCode Extension

<h1 align="center">OpenQC</h1>

<p align="center">
  <strong>Open-Source VSCode Extension for Quantum Chemistry</strong>
</p>

<p align="center">
  <a href="https://code.visualstudio.com/">
    <img src="https://img.shields.io/badge/VSCode-1.85+-blue.svg" alt="VSCode Version">
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License">
  </a>
  <a href="https://github.com/newtontech/OpenQC-VSCode/stargazers">
    <img src="https://img.shields.io/github/stars/newtontech/OpenQC-VSCode.svg" alt="GitHub Stars">
  </a>
</p>

This extension provides comprehensive tools for quantum chemistry workflows in VSCode, including visualization, editing, format conversion, and AI-assisted modification for all major quantum chemistry input file formats.

You can expect from this extension:
- **Multi-format Support**: Gaussian, VASP, Quantum ESPRESSO, ORCA, CP2K, and more
- **3D Visualization**: Interactive molecular structure viewer
- **Smart Editing**: Syntax highlighting, auto-completion, and parameter validation
- **AI Integration**: Native MCP server for Claude Code and Copilot enhancement
- **Remote Computing**: SSH and Slurm integration for high-performance computing

## Prerequisites

- [Node.js](https://nodejs.org/en/) 18+ and [npm](https://nodejs.org/en/)
- [VS Code](https://code.visualstudio.com/) 1.85+
- Python 3.8+ (for core library features)

## Getting Started

- Run `git clone https://github.com/newtontech/OpenQC-VSCode.git`
- Run `cd OpenQC-VSCode/vscode-extension`
- Run `npm install` in the terminal
- Press `F5` to run the extension in debug mode
- Open a quantum chemistry input file to start using OpenQC

## Features

<!-- FEATURES_BEGIN -->
| Feature | Description | Status |
| ------- | ----------- | ------ |
| **Gaussian Support** | `.gjf`, `.com` file editing and visualization | ✅ Full Support |
| **VASP Support** | `POSCAR`, `CONTCAR`, `INCAR` editing | ✅ Full Support |
| **Quantum ESPRESSO** | `.pw`, `.scf`, `.nscf` file support | ✅ Full Support |
| **ORCA Support** | `.inp` file editing | ✅ Full Support |
| **CP2K Support** | `.inp` file editing | ✅ Full Support |
| **Q-Chem Support** | `.qcinp` file editing | 🔄 In Progress |
| **NWChem Support** | `.nw` file editing | 🔄 In Progress |
| **Psi4 Support** | `.dat` file editing | 📋 Planned |
<!-- FEATURES_END -->

## Extension Capabilities

### Visualization
- [3D Molecular Viewer](src/views/molecular-viewer.ts): Interactive visualization using [window.createWebviewPanel](https://code.visualstudio.com/api/references/vscode-api#window.createWebviewPanel)
- [Crystal Structure View](src/views/crystal-viewer.ts): Unit cell and supercell rendering
- [Orbital Visualization](src/views/orbital-viewer.ts): Molecular orbitals and electron density
- [Trajectory Player](src/views/trajectory-player.ts): MD trajectory animation

### Smart Editing
- [Syntax Highlighting](syntaxes/): TextMate grammars for all QC formats
- [Auto-completion](src/providers/completion-provider.ts): [languages.registerCompletionItemProvider](https://code.visualstudio.com/api/references/vscode-api#languages.registerCompletionItemProvider)
- [Parameter Validation](src/providers/validation-provider.ts): Real-time error checking
- [Code Actions](src/providers/code-action-provider.ts): [languages.registerCodeActionsProvider](https://code.visualstudio.com/api/references/vscode-api#languages.registerCodeActionsProvider)

### AI Integration
- **MCP Server**: Full Model Context Protocol implementation for Claude Code
- **Copilot Enhancement**: AI-assisted parameter suggestions
- **Natural Language**: "Change functional to B3LYP" commands
- **Error Analysis**: AI-powered error interpretation

### Remote Computing
- **SSH Handler**([src/remote/ssh-handler.ts](src/remote/ssh-handler.ts)): One-click server connection
- **Slurm Integration**([src/remote/slurm-interface.ts](src/remote/slurm-interface.ts)): Job submission and monitoring
- **File Transfer**: Drag-and-drop remote file operations

## Installation

### From VSCode Marketplace (Coming Soon)
```bash
code --install-extension newtontech.openqc
```

### From Source
```bash
# Clone the repository
git clone https://github.com/newtontech/OpenQC-VSCode.git
cd OpenQC-VSCode

# Install dependencies
cd vscode-extension
npm install

# Build and package
npm run compile
vsce package

# Install in VSCode
code --install-extension openqc-*.vsix
```

## Usage

1. Open a quantum chemistry input file (`.gjf`, `.com`, `POSCAR`, etc.)
2. Press `Ctrl+Shift+P` → "OpenQC: Visualize"
3. Explore your molecular structure in 3D!

### Quick Commands

| Command | Keybinding | Description |
|---------|------------|-------------|
| OpenQC: Visualize | `Ctrl+Shift+V` | Open 3D molecular viewer |
| OpenQC: Convert Format | - | Convert between QC formats |
| OpenQC: Validate | - | Validate input parameters |
| OpenQC: Submit Job | - | Submit to remote server |

## Project Structure

```
OpenQC-VSCode/
├── vscode-extension/           # VSCode extension (TypeScript)
│   ├── src/
│   │   ├── extension.ts        # Main entry point
│   │   ├── providers/          # Language providers
│   │   │   ├── completion-provider.ts
│   │   │   ├── validation-provider.ts
│   │   │   └── code-action-provider.ts
│   │   ├── views/              # Webview panels
│   │   │   ├── molecular-viewer.ts
│   │   │   ├── crystal-viewer.ts
│   │   │   └── orbital-viewer.ts
│   │   ├── commands/           # Command implementations
│   │   └── remote/             # Remote computing
│   │       ├── ssh-handler.ts
│   │       └── slurm-interface.ts
│   ├── syntaxes/               # TextMate grammars
│   ├── language-configuration/
│   └── package.json
│
├── core/                       # Core Python library
│   └── openqc/
│       ├── parsers/            # File format parsers
│       ├── converters/         # Format converters
│       └── visualizers/        # Visualization engine
│
├── ai-protocols/               # AI protocol implementations
│   ├── mcp-server/             # MCP server for Claude
│   └── acp-adapter/            # ACP protocol adapter
│
├── docs/                       # Documentation
│   ├── user-guide/
│   ├── api-reference/
│   └── examples/
│
├── tests/                      # Test suite
│   ├── unit/
│   ├── integration/
│   └── fixtures/
│
└── examples/                   # Example input files
    ├── gaussian/
    ├── vasp/
    ├── qe/
    └── orca/
```

## AI Integration

### Claude Code (MCP)

OpenQC provides a native MCP server for seamless Claude Code integration:

```json
{
  "mcpServers": {
    "openqc": {
      "command": "python",
      "args": ["${workspaceFolder}/ai-protocols/mcp-server/server.py"]
    }
  }
}
```

**Available MCP Tools:**
- `parse_qc_file` - Parse any QC input file
- `convert_format` - Convert between formats
- `visualize_structure` - Generate 3D visualization
- `modify_parameters` - AI-assisted parameter changes
- `validate_input` - Check for errors
- `submit_job` - Submit to remote servers

### Example Claude Code Session

```
User: Parse this Gaussian file and convert it to VASP format
Claude: [Uses openqc.parse_qc_file, then openqc.convert_format]
Here's your VASP POSCAR file. The structure contains 24 atoms...

User: Change the functional to PBE0 and increase the basis set
Claude: [Uses openqc.modify_parameters]
Done! I've updated the functional to PBE0 and changed the basis set to def2-TZVP.
```

## Documentation

- [User Guide](docs/user-guide/README.md) - Getting started and usage instructions
- [API Reference](docs/api-reference/README.md) - Extension API documentation
- [Contributing](CONTRIBUTING.md) - How to contribute to OpenQC
- [AGENTS.md](AGENTS.md) - AI agent configuration and workflows

## API & Contribution Points

This extension contributes the following VS Code APIs:

- [window.createWebviewPanel](https://code.visualstudio.com/api/references/vscode-api#window.createWebviewPanel) - 3D molecular visualization
- [languages.registerCompletionItemProvider](https://code.visualstudio.com/api/references/vscode-api#languages.registerCompletionItemProvider) - Auto-completion
- [languages.registerCodeActionsProvider](https://code.visualstudio.com/api/references/vscode-api#languages.registerCodeActionsProvider) - Code actions
- [commands.registerCommand](https://code.visualstudio.com/api/references/vscode-api#commands.registerCommand) - Extension commands
- [contributes.languages](https://code.visualstudio.com/api/references/contribution-points#contributes.languages) - Language support
- [contributes.views](https://code.visualstudio.com/api/references/contribution-points#contributes.views) - Custom views
- [contributes.menus](https://code.visualstudio.com/api/references/contribution-points#contributes.menus) - Context menus

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

1. Fork the repository
2. Clone your fork: `git clone https://github.com/YOUR_USERNAME/OpenQC-VSCode.git`
3. Install dependencies: `cd vscode-extension && npm install`
4. Open in VS Code and press `F5` to run

### Running Tests

```bash
cd vscode-extension
npm test
```

## Community

- [Discussions](https://github.com/newtontech/OpenQC-VSCode/discussions) - Ask questions and share ideas
- [Issues](https://github.com/newtontech/OpenQC-VSCode/issues) - Report bugs and request features
- [AGENTS.md](AGENTS.md) - AI agent workflows and best practices

## License

[MIT](LICENSE) - Copyright (c) 2024 OpenQC Team

---

<p align="center">
  Made with ❤️ by the OpenQC Team
</p>
