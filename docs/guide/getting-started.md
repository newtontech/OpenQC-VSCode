# Getting Started

OpenQC-VSCode is a universal VS Code extension designed for quantum chemistry researchers and computational chemists. It provides seamless support for multiple quantum chemistry software packages.

## Features

- **Multi-Software Support**: CP2K, VASP, Gaussian, ORCA, Quantum ESPRESSO, GAMESS, NWChem
- **Syntax Highlighting**: Rich, context-aware highlighting for all input file formats
- **IntelliSense**: Smart code completion and parameter suggestions
- **3D Visualization**: Interactive molecular structure visualization
- **Validation**: Real-time error detection and validation
- **LSP Integration**: Full Language Server Protocol support

## Quick Start

### Installation

1. Open VS Code
2. Go to Extensions (Ctrl+Shift+X)
3. Search for "OpenQC-VSCode"
4. Click Install

### Open a Quantum Chemistry File

The extension automatically activates when you open supported file types:

- **CP2K**: `.inp` files
- **VASP**: `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR`
- **Gaussian**: `.gjf`, `.com` files
- **ORCA**: `.inp` files
- **Quantum ESPRESSO**: `.in`, `.pw.in`, `.relax.in`, etc.
- **GAMESS**: `.inp` files
- **NWChem**: `.nw`, `.nwinp` files

### Commands

Open the Command Palette (Ctrl+Shift+P) and type "OpenQC":

- `OpenQC: Visualize Molecular Structure` - Open 3D structure view
- `OpenQC: Plot Calculation Data` - Plot output data
- `OpenQC: Validate Input File` - Validate current file
- `OpenQC: Start Language Server` - Start LSP for current file
- `OpenQC: Restart Language Server` - Restart LSP

## Next Steps

- [Installation](#installation) - Install the extension
- [Commands](#commands) - Use the extension commands
- [Features](#features) - Learn about core features
