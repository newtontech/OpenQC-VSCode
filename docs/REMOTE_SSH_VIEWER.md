# Remote SSH Viewer Workflow

## Overview

OpenQC supports VS Code Remote SSH workflows for viewing molecular and crystal structures on remote machines (HPC clusters, lab servers, etc.) without requiring a GUI session on the remote host.

## Architecture

```text
┌─────────────────────┐     Remote SSH      ┌──────────────────────┐
│   Local Machine      │ ◀─────────────────▶ │   Remote Server      │
│   (VS Code UI)       │                     │   (Extension Host)   │
│                      │                     │                      │
│   Webview (3Dmol.js) │                     │   Python Bridge      │
│   Renders locally    │                     │   Parses files       │
│                      │                     │                      │
│   Uses bundled       │                     │   Reads structures   │
│   extension assets   │                     │   No GUI needed      │
└─────────────────────┘                     └──────────────────────┘
```

### Key Points

1. **Extension host runs on the remote machine**: Commands, Python bridge, and file parsing happen on the remote server.
2. **Webview UI runs locally**: The 3D viewer renders in your local browser using bundled extension assets.
3. **No CDN dependency**: The viewer loads 3Dmol.js from bundled extension files via `webview.asWebviewUri()`.
4. **No internet required**: After extension installation, all viewer assets are local.

## Prerequisites

### Remote Server
- Python 3.8+ (for Python bridge features)
- Optional: `pymatgen`, `ase`, `cclib`, `dpdata` for extended format support

### Local Machine
- VS Code with Remote-SSH extension
- OpenQC extension installed

## Setup

1. Install OpenQC extension in VS Code
2. Connect to remote server via Remote SSH
3. Open a structure file (POSCAR, CIF, XYZ, etc.)
4. Run `OpenQC: Visualize Molecular Structure`

## Python Dependencies

Install on the **remote** machine where the extension host runs:

```bash
# Core (optional but recommended)
pip install numpy

# Structure parsing
pip install pymatgen ase

# Output parsing
pip install cclib

# Dataset support
pip install dpdata
```

Check backend status:
```
OpenQC: Check Scientific Python Backend
```

## Supported Formats

| Format | Native | Python Bridge | Periodic |
|--------|--------|---------------|----------|
| XYZ    | ✅     | ✅            | ❌       |
| POSCAR | ✅     | ✅            | ✅       |
| CIF    | ❌     | ✅ (pymatgen) | ✅       |
| Gaussian input | ✅ | ✅         | ❌       |
| ORCA input     | ✅ | ✅         | ❌       |
| Gaussian output| ❌ | ✅ (cclib)  | ❌       |
| ORCA output    | ❌ | ✅ (cclib)  | ❌       |
| QE input       | ✅ | ✅         | ✅       |
| CP2K input     | ✅ | ✅         | ✅       |

## Troubleshooting

### WebGL Not Available
- Update your local browser/VS Code
- Check GPU acceleration settings in VS Code

### Python Not Found
- Set `openqc.pythonPath` to the full path of Python on the remote server
- Run `OpenQC: Configure Python Path`

### Missing Packages
- Run `OpenQC: Check Scientific Python Backend` to see which packages are available
- Install missing packages on the **remote** machine

### Extension Installed Locally But Not Remotely
- Ensure the extension is installed in the remote VS Code server
- Extensions need to be installed separately for remote contexts

### Large Structures
- Structures with >10,000 atoms may be slow to render
- Supercell previews have a 10,000 atom guardrail
- Use the periodic controls to generate smaller supercells first

## Remote SSH Smoke Test Checklist

- [ ] Connect to remote folder via Remote SSH
- [ ] Open a POSCAR file
- [ ] Run `OpenQC: Visualize Molecular Structure`
- [ ] Verify unit cell renders in the 3D viewer
- [ ] Run `OpenQC: Check Scientific Python Backend`
- [ ] Verify Python version and package status shown
- [ ] Toggle supercell preview (2x2x2)
- [ ] Export a PNG snapshot
- [ ] Test with no internet connection (viewer still works)
