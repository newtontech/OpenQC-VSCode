# AGENTS.md - OpenQC Agent Configuration

This file defines how AI agents should interact with the OpenQC VSCode Extension project.

## Agent Roles

### 1. Code Agent
**Purpose**: Implement features, fix bugs, and maintain code quality

**Responsibilities**:
- Write TypeScript code following the existing patterns
- Implement new quantum chemistry format parsers
- Add new VSCode commands and providers
- Write unit and integration tests

**Best Practices**:
- Follow the existing code style (enforced by ESLint)
- Use the VSCode API consistently with existing code
- Write descriptive commit messages
- Add tests for new features
- Update documentation when adding new features

**Key Files**:
- `vscode-extension/src/extension.ts` - Main entry point
- `vscode-extension/src/providers/` - Language providers
- `vscode-extension/src/views/` - Webview implementations
- `vscode-extension/src/commands/` - Command implementations

**Example Workflow**:
```
1. Read existing similar implementation
2. Create new feature in appropriate directory
3. Register in extension.ts
4. Add tests
5. Update README.md with new feature
6. Create PR with clear description
```

---

### 2. Documentation Agent
**Purpose**: Maintain and improve project documentation

**Responsibilities**:
- Update README.md with new features
- Write user guides and tutorials
- Document API changes
- Maintain AGENTS.md

**Best Practices**:
- Use clear, concise language
- Include code examples
- Update documentation when features change
- Follow VSCode documentation conventions

**Key Files**:
- `README.md` - Main project documentation
- `docs/user-guide/` - User documentation
- `docs/api-reference/` - API documentation
- `CONTRIBUTING.md` - Contribution guidelines

---

### 3. Testing Agent
**Purpose**: Ensure code quality through comprehensive testing

**Responsibilities**:
- Write unit tests for new features
- Create integration tests
- Test across different quantum chemistry formats
- Verify VSCode API compatibility

**Best Practices**:
- Test both success and error cases
- Use realistic test data
- Mock VSCode API appropriately
- Test format conversions thoroughly

**Key Directories**:
- `vscode-extension/src/test/` - Test files
- `tests/unit/` - Unit tests
- `tests/integration/` - Integration tests
- `tests/fixtures/` - Test data files

---

### 4. Integration Agent
**Purpose**: Integrate with external tools and services

**Responsibilities**:
- Implement MCP server features
- Add ACP protocol support
- Integrate with quantum chemistry software
- Maintain remote computing features

**Best Practices**:
- Follow protocol specifications exactly
- Handle errors gracefully
- Document protocol implementations
- Test with real quantum chemistry workflows

**Key Files**:
- `ai-protocols/mcp-server/server.py` - MCP server
- `ai-protocols/acp-adapter/adapter.py` - ACP adapter
- `src/remote/ssh-handler.ts` - SSH connection
- `src/remote/slurm-interface.ts` - Slurm integration

---

## Project Structure for Agents

```
OpenQC-VSCode/
├── vscode-extension/           # Main VSCode extension
│   ├── src/
│   │   ├── extension.ts        # Extension entry point
│   │   ├── providers/          # Language providers (completion, validation)
│   │   ├── views/              # Webview panels (3D visualization)
│   │   ├── commands/           # VSCode commands
│   │   └── remote/             # Remote computing (SSH, Slurm)
│   ├── syntaxes/               # Syntax highlighting definitions
│   ├── src/test/               # Extension tests
│   └── package.json            # Extension manifest
│
├── core/                       # Python core library
│   └── openqc/                 # Main Python package
│       ├── parsers/            # File format parsers
│       ├── converters/         # Format converters
│       └── visualizers/        # 3D visualization engine
│
├── ai-protocols/               # AI protocol implementations
│   ├── mcp-server/             # Model Context Protocol server
│   └── acp-adapter/            # Agent Control Protocol adapter
│
├── docs/                       # Documentation
│   ├── user-guide/             # User documentation
│   └── api-reference/          # API documentation
│
├── tests/                      # Test suite
│   ├── unit/                   # Unit tests
│   ├── integration/            # Integration tests
│   └── fixtures/               # Test data (QC input files)
│
└── examples/                   # Example QC input files
    ├── gaussian/               # Gaussian examples
    ├── vasp/                   # VASP examples
    ├── qe/                     # Quantum ESPRESSO examples
    └── orca/                   # ORCA examples
```

---

## Common Tasks

### Adding Support for New QC Software

1. **Parser** (`core/openqc/parsers/`)
   - Create parser class for the new format
   - Implement `parse()` and `serialize()` methods
   - Add tests with sample files

2. **Syntax Highlighting** (`vscode-extension/syntaxes/`)
   - Create TextMate grammar file
   - Register in `package.json`

3. **Language Provider** (`vscode-extension/src/providers/`)
   - Add completion items
   - Add validation rules
   - Add code actions

4. **Documentation** (`docs/user-guide/`)
   - Add user guide section
   - Update supported formats table in README.md

### Adding New MCP Tools

1. Implement tool in `ai-protocols/mcp-server/server.py`
2. Add tool schema to server capabilities
3. Document tool in AGENTS.md
4. Add example usage in README.md

---

## Code Style

### TypeScript (VSCode Extension)

```typescript
// Use explicit types
function parseInputFile(filePath: string): QCInput {
    // Implementation
}

// Use async/await for async operations
async function visualizeStructure(input: QCInput): Promise<void> {
    await createWebviewPanel(input);
}

// Use interfaces for data structures
interface QCInput {
    software: string;
    atoms: Atom[];
    parameters: Record<string, any>;
}
```

### Python (Core Library)

```python
# Use type hints
def parse_gjf(file_path: str) -> GaussianInput:
    """Parse a Gaussian input file."""
    # Implementation
    pass

# Use dataclasses for data structures
from dataclasses import dataclass

@dataclass
class Atom:
    element: str
    x: float
    y: float
    z: float
```

---

## Testing Guidelines

### Unit Tests

```typescript
// Test individual functions
suite('Parser Tests', () => {
    test('should parse Gaussian input', () => {
        const result = parseGaussianInput(sampleGjf);
        assert.strictEqual(result.software, 'Gaussian');
        assert.strictEqual(result.atoms.length, 24);
    });
});
```

### Integration Tests

```typescript
// Test end-to-end workflows
suite('Extension Integration', () => {
    test('should open and visualize file', async () => {
        const document = await openDocument('test.gjf');
        await executeCommand('openqc.visualize', document);
        // Assert webview is created
    });
});
```

---

## Documentation Standards

### README Updates

When adding features, update:
1. Features table with new capability
2. API & Contribution Points section
3. Usage examples if applicable

### User Guide

Structure:
```markdown
# Feature Name

## Overview
Brief description of the feature.

## Usage
Step-by-step instructions.

## Examples
Code examples and screenshots.

## Troubleshooting
Common issues and solutions.
```

---

## Commit Message Format

```
type(scope): description

[optional body]

[optional footer]
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation only
- `test`: Adding tests
- `refactor`: Code refactoring
- `chore`: Maintenance tasks

Examples:
```
feat(parser): add NWChem input parser

fix(validation): correct parameter validation for ORCA
docs(readme): update supported formats table
test(views): add tests for molecular viewer
```

---

## Resources

- [VSCode Extension API](https://code.visualstudio.com/api/references/vscode-api)
- [VSCode Contribution Points](https://code.visualstudio.com/api/references/contribution-points)
- [MCP Protocol](https://modelcontextprotocol.io/)
- [Quantum Chemistry Software Manuals](docs/user-guide/qc-software.md)

---

Last Updated: 2024-02-28
