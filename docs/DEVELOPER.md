# Developer Guide

This guide is for developers who want to contribute to or extend OpenQC-VSCode.

## Project Structure

```
OpenQC-VSCode/
├── src/                      # TypeScript source code
│   ├── extension.ts          # Extension entry point
│   ├── managers/             # LSP and file management
│   │   ├── LSPManager.ts     # Language server management
│   │   └── FileTypeDetector.ts
│   ├── providers/            # Webview providers
│   │   ├── StructureViewer.ts
│   │   └── DataPlotter.ts
│   └── visualizers/          # Visualization components
│       └── Molecule3D.ts
├── syntaxes/                 # TextMate grammars
├── language-configurations/  # Language configuration files
├── examples/                 # Example input files
├── docs/                     # Documentation
└── tests/                    # Test files
```

## Architecture

### Extension Lifecycle

1. **Activation**: Triggered by opening a supported file type
2. **Detection**: FileTypeDetector identifies the quantum chemistry software
3. **LSP Launch**: LSPManager starts the appropriate language server
4. **Feature Registration**: Commands and providers are registered

### Data Flow

```
File Open
    ↓
FileTypeDetector.detectSoftware()
    ↓
LSPManager.startLSPForDocument()
    ↓
User Command (e.g., Visualize Structure)
    ↓
Molecule3D.parseAtoms()
    ↓
StructureViewer Webview
```

## Adding a New Quantum Chemistry Package

### 1. Add Language Configuration

Create `language-configurations/{package}.json`:

```json
{
  "comments": {
    "lineComment": "!"
  },
  "brackets": [
    ["{", "}"],
    ["[", "]"],
    ["(", ")"]
  ]
}
```

### 2. Add Syntax Grammar

Create `syntaxes/{package}.tmLanguage.json`:

```json
{
  "scopeName": "source.{package}",
  "patterns": [
    {
      "match": "^\\s*#.*$",
      "name": "comment.line.number-sign.{package}"
    }
  ]
}
```

### 3. Update package.json

Add to `contributes.languages`:

```json
{
  "id": "{package}",
  "aliases": ["PackageName", "{package}"],
  "extensions": [".ext"],
  "configuration": "./language-configurations/{package}.json"
}
```

Add to `contributes.grammars`:

```json
{
  "language": "{package}",
  "scopeName": "source.{package}",
  "path": "./syntaxes/{package}.tmLanguage.json"
}
```

### 4. Update FileTypeDetector

Add detection pattern in `src/managers/FileTypeDetector.ts`:

```typescript
{
  software: 'NewPackage',
  extensions: ['.ext'],
  contentPatterns: [
    /specific pattern/i
  ]
}
```

### 5. Update LSPManager

Add to `getLanguageId()` and `getExtensions()` methods.

### 6. Add Parser

Update `src/visualizers/Molecule3D.ts`:

```typescript
private parseNewPackageAtoms(content: string): Atom[] {
  // Implementation
}
```

### 7. Add Tests

Create test files in `tests/fixtures/{package}/`.

## Testing

### Unit Tests

```bash
npm run test:unit
```

### Integration Tests

```bash
npm run test:integration
```

### E2E Tests

```bash
npm run test:e2e
```

## Debugging

1. Open the project in VSCode
2. Press `F5` to launch Extension Development Host
3. Open a quantum chemistry file to test

## Code Style

- Use TypeScript strict mode
- Follow ESLint configuration
- Add JSDoc comments for public APIs
- Write tests for new features

## Submitting Changes

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Update documentation
6. Submit a pull request

## Release Process

1. Update version in `package.json`
2. Update `CHANGELOG.md`
3. Run tests: `npm test`
4. Package: `npx vsce package`
5. Publish: `npx vsce publish`

## Resources

- [VSCode Extension API](https://code.visualstudio.com/api)
- [Language Server Protocol](https://microsoft.github.io/language-server-protocol/)
- [3Dmol.js Documentation](https://3dmol.csb.pitt.edu/doc/index.html)
- [Plotly.js Documentation](https://plotly.com/javascript/)
