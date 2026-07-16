# OpenQC-VSCode API Documentation

## Extension API

### Commands

The extension provides the following commands:

#### `openqc.visualizeStructure`
Opens the 3D molecular structure viewer for the current file.

**Usage:**
```typescript
vscode.commands.executeCommand('openqc.visualizeStructure');
```

#### `openqc.plotData`
Opens the data plotter for the current file.

**Usage:**
```typescript
vscode.commands.executeCommand('openqc.plotData');
```

#### `openqc.previewInput`
Opens the structured input preview for the current file.

**Usage:**
```typescript
vscode.commands.executeCommand('openqc.previewInput');
```

#### `openqc.startLSP`
Manually starts the language server for the current file type.

**Usage:**
```typescript
vscode.commands.executeCommand('openqc.startLSP');
```

#### `openqc.stopLSP`
Stops the language server for the current file type.

**Usage:**
```typescript
vscode.commands.executeCommand('openqc.stopLSP');
```

#### `openqc.restartLSP`
Restarts the language server for the current file type.

**Usage:**
```typescript
vscode.commands.executeCommand('openqc.restartLSP');
```

---

## Configuration API

### LSP Settings

```json
{
  "openqc.lsp.abacus.enabled": true,
  "openqc.lsp.abacus.command": "abacus-lsp",
  "openqc.lsp.abinit.enabled": true,
  "openqc.lsp.abinit.command": "abinit-lsp",
  "openqc.lsp.cif.enabled": true,
  "openqc.lsp.cif.command": "cif-lsp",
  "openqc.lsp.cp2k.enabled": true,
  "openqc.lsp.cp2k.command": "cp2k-language-server",
  "openqc.lsp.vasp.enabled": true,
  "openqc.lsp.vasp.command": "vasp-lsp",
  "openqc.lsp.gaussian.enabled": true,
  "openqc.lsp.gaussian.command": "gaussian-lsp",
  "openqc.lsp.orca.enabled": true,
  "openqc.lsp.orca.command": "orca-lsp",
  "openqc.lsp.qe.enabled": true,
  "openqc.lsp.qe.command": "qe-lsp",
  "openqc.lsp.gamess.enabled": true,
  "openqc.lsp.gamess.command": "gamess-lsp",
  "openqc.lsp.nwchem.enabled": true,
  "openqc.lsp.nwchem.command": "nwchem-lsp",
  "openqc.lsp.gpumd.enabled": true,
  "openqc.lsp.gpumd.command": "gpumd-lsp",
  "openqc.lsp.gromacs.enabled": true,
  "openqc.lsp.gromacs.command": "gromacs-lsp",
  "openqc.lsp.lammps.enabled": true,
  "openqc.lsp.lammps.command": "lmp-lsp",
  "openqc.lsp.mlip.enabled": true,
  "openqc.lsp.mlip.command": "mlip-lsp",
  "openqc.lsp.pyatb.enabled": true,
  "openqc.lsp.pyatb.command": "pyatb-lsp",
  "openqc.lsp.pyscf.enabled": true,
  "openqc.lsp.pyscf.command": "pyscf-lsp"
}
```

`openqc.lsp.<languageId>.path` remains a deprecated alias for older settings. Prefer `openqc.lsp.<languageId>.command` for new configuration.

### Visualization Settings

```json
{
  "openqc.visualization.moleculeRenderer": "3Dmol.js",
  "openqc.visualization.plotLibrary": "Plotly.js",
  "openqc.visualization.autoOpen": true
}
```

---

## Language IDs

| Package | Language ID | File Extensions |
|---------|-------------|-----------------|
| ABACUS | `abacus` | `INPUT`, `STRU`, `KPT` |
| ABINIT | `abinit` | `.abi`, `.abinit` |
| CIF | `cif` | `.cif` |
| CP2K | `cp2k` | `.inp` |
| VASP | `vasp` | `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR` |
| Gaussian | `gaussian` | `.gjf`, `.com` |
| ORCA | `orca` | `.inp` |
| Quantum ESPRESSO | `qe` | `.in`, `.pw.in`, `.relax.in`, etc. |
| GAMESS | `gamess` | `.inp` |
| NWChem | `nwchem` | `.nw`, `.nwinp` |
| GPUMD | `gpumd` | `run.in`, `nep.in` |
| GROMACS | `gromacs` | `.top`, `.itp`, `.mdp`, `.gro` |
| LAMMPS | `lammps` | `.lmp`, `.lammps`, `.lmps`, `in.lammps` |
| MLIP | `mlip` | `.mlip.json`, `.mlip.yaml`, `.mlip.yml` |
| PyATB | `pyatb` | `.pyatb.py`, `run_pyatb.py` |
| PySCF | `pyscf` | `.pyscf.py`, `run_pyscf.py` |

---

## Extension Development

### Building from Source

```bash
git clone https://github.com/newtontech/OpenQC-VSCode.git
cd OpenQC-VSCode
nvm use
npm ci
npm run compile
```

### Running Tests

```bash
npm test
```

### Packaging

```bash
npm run check:release
```

---

## Internal APIs

### LSPManager

Manages language server connections.

```typescript
class LSPManager {
  startLSPForDocument(document: vscode.TextDocument): Promise<void>;
  stopLSPForDocument(document: vscode.TextDocument): Promise<void>;
  restartLSPForDocument(document: vscode.TextDocument): Promise<void>;
}
```

### FileTypeDetector

Detects quantum chemistry software from file content.

```typescript
class FileTypeDetector {
  detectSoftware(document: vscode.TextDocument): QuantumChemistrySoftware | null;
  getSoftwareInfo(software: QuantumChemistrySoftware): SoftwareInfo;
}
```

### Molecule3D

Parses atomic coordinates from various formats.

```typescript
class Molecule3D {
  parseAtoms(content: string, software: QuantumChemistrySoftware): Atom[];
}
```

---

## Contributing

See [CONTRIBUTING.md](https://github.com/newtontech/OpenQC-VSCode/blob/master/CONTRIBUTING.md) for development guidelines.
