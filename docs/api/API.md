# API Reference

This document provides reference documentation for the OpenQC-VSCode API.

## Table of Contents

- [Extension API](#extension-api)
- [Parser API](#parser-api)
- [Converter API](#converter-api)
- [Visualization API](#visualization-api)
- [AI API](#ai-api)

---

## Extension API

### Commands

All commands are accessible via VSCode command palette (`Ctrl+Shift+P` / `Cmd+Shift+P`).

| Command | Description |
|---------|-------------|
| `openqc.parse` | Parse current file |
| `openqc.visualize` | Open 3D visualization |
| `openqc.convert` | Convert file to different format |
| `openqc.validate` | Validate file for errors |
| `openqc.ai.suggest` | Get AI suggestions |
| `openqc.ai.explain` | Explain parameters |

### Configuration

Settings in `settings.json`:

```json
{
  "openqc.ai.provider": "openai",
  "openqc.ai.model": "gpt-4",
  "openqc.visualization.renderer": "ngl",
  "openqc.parsers.preferPython": false
}
```

---

## Parser API

### Base Parser Interface

```typescript
interface Parser {
  /**
   * Parse file content into structured data
   * @param content - Raw file content
   * @returns Parsed result with data and metadata
   */
  parse(content: string): ParseResult;
  
  /**
   * Validate parsed result
   * @param result - Parse result to validate
   * @returns Validation result with errors/warnings
   */
  validate(result: ParseResult): ValidationResult;
  
  /**
   * Get syntax highlighting tokens
   * @param content - Raw file content
   * @returns Array of tokens with positions
   */
  tokenize(content: string): Token[];
}

interface ParseResult {
  success: boolean;
  data: any;
  metadata: {
    format: string;
    version?: string;
    warnings?: string[];
  };
  errors?: ParseError[];
}

interface ValidationResult {
  valid: boolean;
  errors: ValidationError[];
  warnings: ValidationWarning[];
}
```

### INCAR Parser

```typescript
class IncarParser implements Parser {
  /**
   * Parse VASP INCAR file
   * @example
   * const parser = new IncarParser();
   * const result = parser.parse('ENCUT = 520\nISMEAR = 0');
   * // result.data = { ENCUT: 520, ISMEAR: 0 }
   */
  parse(content: string): ParseResult;
  
  /**
   * Convert parsed INCAR back to string
   */
  stringify(data: Record<string, any>): string;
  
  /**
   * Get documentation for a parameter
   */
  getParameterDocs(param: string): string | undefined;
}
```

### POSCAR Parser

```typescript
class PoscarParser implements Parser {
  /**
   * Parse VASP POSCAR/CONTCAR file
   * @example
   * const parser = new PoscarParser();
   * const result = parser.parse(poscarContent);
   * // result.data = { lattice, atoms, positions, ... }
   */
  parse(content: string): ParseResult;
  
  /**
   * Get atomic structure in standard format
   */
  getStructure(result: ParseResult): AtomicStructure;
}

interface AtomicStructure {
  lattice: Matrix3x3;
  atoms: Atom[];
  positions: Vector3[];
  species: string[];
  cellType?: 'direct' | 'cartesian';
}

interface Atom {
  element: string;
  position: Vector3;
  mass?: number;
  charge?: number;
}
```

### Gaussian Parser

```typescript
class GaussianParser implements Parser {
  /**
   * Parse Gaussian input file (.com, .gjf)
   */
  parse(content: string): ParseResult;
  
  /**
   * Get route section
   */
  getRouteSection(result: ParseResult): string;
  
  /**
   * Get molecular specification
   */
  getMolecule(result: ParseResult): Molecule;
}

interface Molecule {
  charge: number;
  multiplicity: number;
  atoms: Atom[];
}
```

---

## Converter API

### Base Converter Interface

```typescript
interface Converter {
  /**
   * Convert from one format to another
   * @param input - Source content
   * @param from - Source format
   * @param to - Target format
   * @returns Converted content
   */
  convert(input: string, from: Format, to: Format): Promise<string>;
  
  /**
   * Check if conversion is supported
   */
  canConvert(from: Format, to: Format): boolean;
  
  /**
   * Get list of supported target formats
   */
  getSupportedTargets(from: Format): Format[];
}

type Format = 
  | 'vasp-poscar'
  | 'vasp-incar'
  | 'gaussian'
  | 'orca'
  | 'xyz'
  | 'pdb'
  | 'cif';
```

### VASP Converter

```typescript
class VASPConverter implements Converter {
  /**
   * Convert VASP POSCAR to other formats
   * @example
   * const converter = new VASPConverter();
   * const gaussian = await converter.convert(
   *   poscarContent,
   *   'vasp-poscar',
   *   'gaussian'
   * );
   */
  convert(input: string, from: Format, to: Format): Promise<string>;
  
  /**
   * Batch convert multiple files
   */
  convertBatch(files: File[]): Promise<Map<string, string>>;
}
```

### Universal Converter

```typescript
class UniversalConverter {
  /**
   * Auto-detect format and convert
   */
  autoConvert(input: string, to: Format): Promise<string>;
  
  /**
   * Get conversion chain for unsupported direct conversions
   */
  getConversionPath(from: Format, to: Format): Format[];
}
```

---

## Visualization API

### 3D Viewer

```typescript
interface Viewer3D {
  /**
   * Load molecular structure
   */
  loadStructure(structure: AtomicStructure): void;
  
  /**
   * Set rendering mode
   */
  setRepresentation(mode: RepresentationMode): void;
  
  /**
   * Highlight specific atoms
   */
  highlightAtoms(indices: number[]): void;
  
  /**
   * Get camera state
   */
  getCameraState(): CameraState;
  
  /**
   * Set camera position
   */
  setCamera(state: CameraState): void;
}

type RepresentationMode = 
  | 'ball-stick'
  | 'space-filling'
  | 'wireframe'
  | 'cartoon';
```

### Webview Messaging

```typescript
// Extension → Webview
interface ExtensionMessage {
  type: 'load' | 'update' | 'highlight';
  payload: any;
}

// Webview → Extension
interface WebviewMessage {
  type: 'ready' | 'selection' | 'edit';
  payload: any;
}

class WebviewMessaging {
  /**
   * Send message to webview
   */
  postMessage(message: ExtensionMessage): void;
  
  /**
   * Listen for webview messages
   */
  onMessage(callback: (msg: WebviewMessage) => void): Disposable;
}
```

---

## AI API

### AI Provider Interface

```typescript
interface AIProvider {
  /**
   * Generate completion
   */
  complete(prompt: string, options?: CompletionOptions): Promise<string>;
  
  /**
   * Stream completion
   */
  streamComplete(prompt: string, onChunk: (chunk: string) => void): Promise<void>;
}

interface CompletionOptions {
  model?: string;
  temperature?: number;
  maxTokens?: number;
  context?: string[];
}
```

### AI Service

```typescript
class AIService {
  /**
   * Analyze input file and suggest improvements
   */
  analyzeFile(content: string, format: Format): Promise<AnalysisResult>;
  
  /**
   * Explain parameters in plain English
   */
  explainParameters(params: Record<string, any>): Promise<string>;
  
  /**
   * Generate input file from description
   */
  generateFromDescription(description: string, format: Format): Promise<string>;
  
  /**
   * Debug failed calculation
   */
  debugError(error: string, context: FileContext): Promise<DebugSuggestion>;
}

interface AnalysisResult {
  suggestions: Suggestion[];
  warnings: string[];
  optimizations: Optimization[];
}

interface Suggestion {
  parameter: string;
  currentValue: any;
  suggestedValue: any;
  reason: string;
  confidence: number;
}
```

---

## Error Handling

### Custom Errors

```typescript
class ParseError extends Error {
  constructor(
    message: string,
    public line?: number,
    public column?: number
  ) {
    super(message);
  }
}

class ConversionError extends Error {
  constructor(
    message: string,
    public from: Format,
    public to: Format
  ) {
    super(message);
  }
}

class AIError extends Error {
  constructor(
    message: string,
    public provider: string
  ) {
    super(message);
  }
}
```

---

## Events

### Extension Events

```typescript
// File parsed event
vscode.commands.registerCommand('openqc.onFileParsed', (result: ParseResult) => {
  // Handle parsed file
});

// Structure updated event
vscode.commands.registerCommand('openqc.onStructureUpdated', (structure: AtomicStructure) => {
  // Handle structure update
});
```

---

## TypeScript Types

All TypeScript type definitions are exported from `@openqc/types`:

```typescript
import type {
  Parser,
  ParseResult,
  Converter,
  AtomicStructure,
  Viewer3D,
  AIProvider
} from '@openqc/types';
```

---

## Python API

For Python integrations, see the Python API documentation:

```python
from openqc.adapters import DPDataAdapter

# Convert file
adapter = DPDataAdapter()
result = adapter.convert(poscar_content, 'vasp', 'gaussian')
```

---

## Versioning

The API follows semantic versioning:
- **Major**: Breaking changes
- **Minor**: New features, backward compatible
- **Patch**: Bug fixes

Current API version: `1.0.0`

---

## Future API

Planned features for future versions:

- WebSocket API for real-time updates
- REST API for external tools
- Plugin API for custom extensions
- Batch processing API

---

**Last Updated**: 2026-02-28
**Version**: 1.0.0
