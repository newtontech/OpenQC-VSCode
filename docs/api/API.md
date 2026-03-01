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
| CP2K | `cp2k` | `.inp` |
| VASP | `vasp` | `INCAR`, `POSCAR`, `KPOINTS`, `POTCAR` |
| Gaussian | `gaussian` | `.gjf`, `.com` |
| ORCA | `orca` | `.inp` |
| Quantum ESPRESSO | `qe` | `.in`, `.pw.in`, `.relax.in`, etc. |
| GAMESS | `gamess` | `.inp` |
| NWChem | `nwchem` | `.nw`, `.nwinp` |

---

## Extension Development

### Building from Source

```bash
git clone https://github.com/newtontech/OpenQC-VSCode.git
cd OpenQC-VSCode
npm install
npm run compile
```

### Running Tests

```bash
npm test
```

### Packaging

```bash
npx vsce package
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

See [CONTRIBUTING.md](../CONTRIBUTING.md) for development guidelines.
