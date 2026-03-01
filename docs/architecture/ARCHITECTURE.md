# Architecture Overview

This document describes the high-level architecture of OpenQC-VSCode.

## System Architecture

OpenQC-VSCode follows a layered architecture with clear separation of concerns:

```
┌─────────────────────────────────────────────────────────┐
│                    VSCode Extension                      │
│  ┌──────────────────────────────────────────────────┐  │
│  │              UI Layer (React Webview)             │  │
│  │  • 3D Visualization (Three.js/NGL)                │  │
│  │  • Input Forms                                     │  │
│  │  • Sidebar Panels                                  │  │
│  └──────────────────────────────────────────────────┘  │
│  ┌──────────────────────────────────────────────────┐  │
│  │           Business Logic Layer                    │  │
│  │  • Commands                                        │  │
│  │  • Controllers                                     │  │
│  │  • State Management (Zustand/Redux)               │  │
│  └──────────────────────────────────────────────────┘  │
│  ┌──────────────────────────────────────────────────┐  │
│  │              Service Layer                        │  │
│  │  • Parser Service                                  │  │
│  │  • Converter Service                               │  │
│  │  • AI Service                                      │  │
│  └──────────────────────────────────────────────────┘  │
│  ┌──────────────────────────────────────────────────┐  │
│  │               Core Layer                          │  │
│  │  • Parsers (Chevrotain)                            │  │
│  │  • Converters                                      │  │
│  │  • Validators                                      │  │
│  └──────────────────────────────────────────────────┘  │
│  ┌──────────────────────────────────────────────────┐  │
│  │         External Integrations Layer               │  │
│  │  • dpdata (Python)                                 │  │
│  │  • ASE (Python)                                    │  │
│  │  • PyMOL (Python)                                  │  │
│  │  • AI APIs (OpenAI/Ollama)                         │  │
│  └──────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────┘
```

## Core Components

### 1. Parsers

**Purpose**: Parse quantum chemistry input files into structured data.

**Implementation**:
- TypeScript parsers using Chevrotain (fast, incremental parsing)
- Python fallback using Pyparsing (for complex formats)
- Abstract parser interface for extensibility

**Key Classes**:
```
src/parsers/
├── base/
│   ├── Parser.ts              # Abstract base parser
│   └── ParseResult.ts         # Common result types
├── vasp/
│   ├── IncarParser.ts
│   ├── PoscarParser.ts
│   ├── KpointsParser.ts
│   └── PotcarParser.ts
├── gaussian/
│   └── GaussianParser.ts
├── orca/
│   └── OrcaParser.ts
└── index.ts                   # Parser factory
```

### 2. Converters

**Purpose**: Convert between different file formats.

**Implementation**:
- Primary: dpdata (Python backend)
- TypeScript wrappers for common conversions
- Batch conversion support

**Key Classes**:
```
src/converters/
├── Converter.ts               # Abstract converter
├── VASPConverter.ts
├── GaussianConverter.ts
├── ORCAConverter.ts
├── UniversalConverter.ts      # Multi-format converter
└── adapters/
    ├── DPDataAdapter.ts       # Python integration
    └── ASEAdapter.ts
```

### 3. Visualization

**Purpose**: Render 3D molecular structures.

**Implementation**:
- NGL Viewer or Three.js for rendering
- VSCode Webview for UI
- Real-time updates via bidirectional messaging

**Key Components**:
```
src/visualization/
├── renderer/
│   ├── MoleculeRenderer.ts
│   ├── UnitCellRenderer.ts
│   └── CameraController.ts
├── webview/
│   ├── index.tsx              # React entry point
│   ├── components/
│   │   ├── Viewer3D.tsx
│   │   ├── Toolbar.tsx
│   │   └── InfoPanel.tsx
│   └── state/
│       └── viewerStore.ts     # Zustand store
└── messaging/
    └── WebviewMessaging.ts    # Extension ↔ Webview
```

### 4. AI Service

**Purpose**: Provide AI-powered features.

**Implementation**:
- API clients for OpenAI, Anthropic, Ollama
- Context management (RAG with documentation)
- Prompt engineering for quantum chemistry

**Key Components**:
```
src/ai/
├── providers/
│   ├── OpenAIProvider.ts
│   ├── AnthropicProvider.ts
│   └── OllamaProvider.ts
├── context/
│   ├── ContextBuilder.ts
│   └── RAGSystem.ts
├── prompts/
│   └── quantum-prompts.ts
└── AIService.ts
```

## Data Flow

### Parsing Workflow

```
User opens file
    ↓
VSCode detects file type
    ↓
Parser Factory creates appropriate parser
    ↓
Parser.parse() → ParseResult
    ↓
    ├→ Syntax Highlighting
    ├→ Validation (errors/warnings)
    └→ Data extraction (structure, parameters)
```

### Conversion Workflow

```
User requests conversion
    ↓
Converter validates input
    ↓
Python backend (dpdata) performs conversion
    ↓
Result returned to TypeScript
    ↓
    ├→ Display in editor
    └→ Save to file
```

### Visualization Workflow

```
User opens 3D view
    ↓
Parse structure from file
    ↓
Send to Webview via messaging
    ↓
Webview renders with NGL/Three.js
    ↓
User interacts (rotate, select, edit)
    ↓
Changes sent back to extension
    ↓
Update source file
```

## Design Patterns

### 1. Strategy Pattern (Parsers)

Different parsing strategies for different formats:

```typescript
interface Parser {
  parse(content: string): ParseResult;
  validate(result: ParseResult): ValidationResult;
}

class IncarParser implements Parser { /* ... */ }
class PoscarParser implements Parser { /* ... */ }
```

### 2. Factory Pattern (Parser Creation)

```typescript
class ParserFactory {
  static create(format: FileFormat): Parser {
    switch(format) {
      case 'vasp-incar': return new IncarParser();
      case 'vasp-poscar': return new PoscarParser();
      // ...
    }
  }
}
```

### 3. Adapter Pattern (External Tools)

```typescript
interface ConversionAdapter {
  convert(input: string, from: Format, to: Format): Promise<string>;
}

class DPDataAdapter implements ConversionAdapter {
  async convert(input, from, to) {
    // Call Python dpdata via subprocess
  }
}
```

### 4. Observer Pattern (State Management)

```typescript
// Zustand store for reactive updates
const useViewerStore = create((set) => ({
  structure: null,
  setStructure: (structure) => set({ structure }),
}));
```

## Extension Points

### Adding a New Parser

1. Create parser class implementing `Parser` interface
2. Add to `ParserFactory`
3. Register file type in `package.json`
4. Add syntax highlighting grammar
5. Write tests

### Adding a New Converter

1. Create converter class implementing `Converter` interface
2. Add to `UniversalConverter`
3. Implement Python adapter (if using dpdata)
4. Add UI commands
5. Write tests

### Adding a New AI Provider

1. Create provider class implementing `AIProvider` interface
2. Add configuration in settings
3. Implement in `AIService`
4. Write tests with mocking

## Performance Considerations

### Large Files

- **Lazy loading**: Only parse visible portions
- **WebWorkers**: Heavy computation off main thread
- **Incremental parsing**: Re-parse only changed sections

### 3D Rendering

- **Level of Detail (LOD)**: Reduce detail for distant atoms
- **Instanced rendering**: Efficient for repeated elements
- **Frustum culling**: Only render visible atoms

### Memory Management

- **Dispose patterns**: Clean up Three.js objects
- **Weak references**: For cached structures
- **Streaming**: For very large files

## Security

### External Process Execution

- Validate all inputs before passing to Python
- Sandbox Python execution
- Timeout long-running operations

### AI API Keys

- Store in VSCode secret storage
- Never log or expose keys
- Support local models (Ollama) for air-gapped environments

### File System Access

- Only access files user explicitly opens
- Validate file paths
- Restrict to workspace or explicit directories

## Testing Strategy

### Unit Tests
- Test each parser in isolation
- Mock external dependencies
- Test edge cases and errors

### Integration Tests
- Test Python adapters
- Test VSCode API interactions
- Test file system operations

### E2E Tests
- Test complete workflows
- Test in real VSCode instance
- Cross-platform testing

## Future Considerations

### Extensibility
- Plugin API for custom formats
- Community extensions
- Custom workflows

### Scalability
- Support for very large systems (>10k atoms)
- Multi-file projects
- Distributed computation

### Collaboration
- Real-time collaborative editing
- Shared visualization sessions
- Cloud-based rendering

---

**Last Updated**: 2026-02-28
**Version**: 1.0.0
