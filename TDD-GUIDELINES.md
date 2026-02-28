# Test-Driven Development (TDD) Guidelines

> Best practices for writing tests in OpenQC-VSCode project

## Philosophy

We follow Test-Driven Development (TDD) to ensure code quality, maintainability, and confidence in refactoring. Tests are not just for catching bugs—they're a design tool that helps us write better code.

### Core Principles

1. **Red → Green → Refactor**
   - Write a failing test first
   - Write the minimum code to make it pass
   - Refactor while keeping tests green

2. **Test Behavior, Not Implementation**
   - Focus on what the code does, not how it does it
   - Tests should survive refactoring

3. **Fast Feedback Loop**
   - Tests must run quickly
   - Developers should run tests frequently
   - CI runs all tests on every push

---

## Testing Stack

### TypeScript/VSCode Extension

- **Framework**: Jest
- **Utilities**: @vscode/test-electron
- **Coverage**: Istanbul/NYC
- **Mocking**: jest.mock() + manual mocks

### Python Backend

- **Framework**: pytest
- **Fixtures**: pytest fixtures
- **Coverage**: pytest-cov
- **Mocking**: unittest.mock / pytest-mock

---

## Test Categories

### 1. Unit Tests

**Purpose**: Test individual functions, classes, or modules in isolation.

**Location**: 
- TypeScript: `tests/unit/**/*.test.ts`
- Python: `tests/unit/test_*.py`

**Guidelines**:
- ✅ Fast execution (< 100ms per test)
- ✅ No external dependencies (mock everything)
- ✅ Single responsibility (one assertion per test)
- ✅ Descriptive names: `should_parse_incar_file_correctly`

**Example (TypeScript)**:
```typescript
// tests/unit/parsers/vasp/incarParser.test.ts
import { IncarParser } from '@/parsers/vasp';

describe('IncarParser', () => {
  describe('parse()', () => {
    it('should parse basic INCAR parameters', () => {
      const input = 'ENCUT = 520\nISMEAR = 0';
      const result = IncarParser.parse(input);
      
      expect(result).toEqual({
        ENCUT: 520,
        ISMEAR: 0
      });
    });

    it('should throw error on invalid format', () => {
      const input = 'INVALID_PARAM = abc';
      expect(() => IncarParser.parse(input)).toThrow(ParseError);
    });
  });
});
```

**Example (Python)**:
```python
# tests/unit/test_dpdata_adapter.py
import pytest
from openqc.adapters import DPDataAdapter

def test_convert_vasp_to_gaussian():
    """Should convert VASP POSCAR to Gaussian input format."""
    adapter = DPDataAdapter()
    poscar_content = """
    Test structure
    1.0
    4.0 0.0 0.0
    0.0 4.0 0.0
    0.0 0.0 4.0
    Si
    1
    Direct
    0.0 0.0 0.0
    """
    
    result = adapter.convert(poscar_content, 'vasp', 'gaussian')
    
    assert result.startswith('#')
    assert 'Si' in result
```

---

### 2. Integration Tests

**Purpose**: Test interactions between components or with external systems.

**Location**:
- TypeScript: `tests/integration/**/*.test.ts`
- Python: `tests/integration/test_*.py`

**Guidelines**:
- ✅ Test real interactions (database, file system, APIs)
- ✅ Use test fixtures and sample data
- ✅ Clean up after tests
- ✅ Can be slower (seconds, not milliseconds)

**Example**:
```typescript
// tests/integration/extension.test.ts
import * as vscode from 'vscode';
import { expect } from 'chai';

describe('Extension Integration', () => {
  it('should activate and register commands', async () => {
    const ext = vscode.extensions.getExtension('newtontech.openqc-vscode');
    await ext.activate();
    
    const commands = await vscode.commands.getCommands(true);
    expect(commands).to.include('openqc.parseINCAR');
    expect(commands).to.include('openqc.visualize');
  });
});
```

---

### 3. End-to-End (E2E) Tests

**Purpose**: Test complete user workflows through VSCode UI.

**Location**: `tests/e2e/**/*.test.ts`

**Guidelines**:
- ✅ Test real user scenarios
- ✅ Use VSCode testing framework
- ✅ Slow but comprehensive
- ✅ Run in CI only (not locally)

**Example**:
```typescript
// tests/e2e/workflows.test.ts
import { VSCCode } from '@vscode/test-electron';

describe('User Workflows', () => {
  it('should parse and visualize POSCAR file', async () => {
    // Open file
    await vscode.workspace.openTextDocument('samples/POSCAR');
    
    // Execute command
    await vscode.commands.executeCommand('openqc.visualize');
    
    // Verify webview opens
    const webviews = vscode.window.webviews;
    expect(webviews).to.have.length(1);
  });
});
```

---

## Test-Driven Workflow

### Step-by-Step Process

1. **Write Test First**
   ```bash
   # Create test file
   touch tests/unit/parsers/vasp/poscarParser.test.ts
   
   # Write failing test
   # Run: npm test
   # Expected: FAIL (test doesn't exist yet)
   ```

2. **Implement Code**
   ```bash
   # Create implementation
   touch src/parsers/vasp/poscarParser.ts
   
   # Write minimal code to pass
   # Run: npm test
   # Expected: PASS
   ```

3. **Refactor**
   ```bash
   # Improve code quality
   # Run: npm test (continuously)
   # Expected: Still PASS
   ```

4. **Commit**
   ```bash
   git add .
   git commit -m "feat: add POSCAR parser with tests"
   ```

---

## Coverage Requirements

### Minimum Coverage
- **Overall**: 80%
- **Critical Paths** (parsers, converters): 90%
- **UI Components**: 70%

### Coverage Commands
```bash
# TypeScript
npm run test:coverage

# Python
pytest --cov=src --cov-report=html
```

### Coverage Reports
- Generated in `coverage/` directory
- HTML report: `coverage/lcov-report/index.html`
- CI fails if coverage drops below threshold

---

## Mocking Strategies

### VSCode API Mocking

```typescript
// tests/__mocks__/vscode.ts
export const vscode = {
  window: {
    showInformationMessage: jest.fn(),
    createWebviewPanel: jest.fn(),
  },
  workspace: {
    openTextDocument: jest.fn(),
    getConfiguration: jest.fn(),
  },
  commands: {
    registerCommand: jest.fn(),
    executeCommand: jest.fn(),
  },
};
```

### Python External Tools

```python
# tests/conftest.py
import pytest
from unittest.mock import Mock, patch

@pytest.fixture
def mock_subprocess():
    with patch('subprocess.run') as mock_run:
        mock_run.return_value = Mock(returncode=0, stdout='OK')
        yield mock_run
```

---

## Test Data Management

### Sample Files
- Store in `tests/fixtures/`
- Organize by format: `tests/fixtures/vasp/`, `tests/fixtures/gaussian/`
- Use real-world examples (anonymized if needed)

### Fixtures
```typescript
// tests/fixtures/index.ts
export const SAMPLE_POSCAR = `
Si diamond
5.43
0.5 0.5 0.0
0.0 0.5 0.5
0.5 0.0 0.5
Si
2
Direct
0.00 0.00 0.00
0.25 0.25 0.25
`;

export const INVALID_POSCAR = `
This is not a valid POSCAR
Random text here
`;
```

---

## CI/CD Integration

### GitHub Actions Workflow

```yaml
# .github/workflows/test.yml
name: Tests

on: [push, pull_request]

jobs:
  test-typescript:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-node@v3
        with:
          node-version: 18
      - run: npm ci
      - run: npm run lint
      - run: npm run test:unit
      - run: npm run test:coverage
      - uses: codecov/codecov-action@v3

  test-python:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: 3.11
      - run: pip install -e .[dev]
      - run: pytest
      - run: pytest --cov=src --cov-fail-under=80
```

---

## Best Practices

### Do's ✅

- **Write descriptive test names**
  ```typescript
  it('should throw ParseError when ENCUT is negative', () => {
    // ...
  });
  ```

- **Test edge cases**
  - Empty files
  - Large files
  - Malformed input
  - Unicode characters

- **Use parameterized tests**
  ```python
  @pytest.mark.parametrize("format,expected", [
    ("vasp", "POSCAR"),
    ("gaussian", "input.com"),
    ("orca", "input.inp"),
  ])
  def test_file_extensions(format, expected):
      assert get_default_filename(format) == expected
  ```

- **Clean up resources**
  ```typescript
  afterEach(() => {
    vscode.window.activeTextEditor?.hide();
  });
  ```

### Don'ts ❌

- **Don't test implementation details**
  ```typescript
  // ❌ Bad
  expect(parser.internalState).toBe('ready');
  
  // ✅ Good
  expect(parser.parse('input')).toBeDefined();
  ```

- **Don't use magic numbers**
  ```typescript
  // ❌ Bad
  expect(result.length).toBe(3);
  
  // ✅ Good
  const EXPECTED_ATOM_COUNT = 3;
  expect(result.length).toBe(EXPECTED_ATOM_COUNT);
  ```

- **Don't ignore flaky tests**
  - Fix them immediately
  - Or mark as `.skip()` with a TODO

---

## Testing Checklist

Before submitting a PR:

- [ ] All new code has tests
- [ ] All tests pass locally
- [ ] Coverage ≥ 80%
- [ ] No skipped tests without justification
- [ ] Test names are descriptive
- [ ] Edge cases covered
- [ ] CI passes on all platforms

---

## Tools & Commands

### TypeScript

```bash
# Run all tests
npm test

# Run specific test file
npm test -- incarParser.test.ts

# Run with coverage
npm run test:coverage

# Watch mode
npm run test:watch

# Debug tests
npm run test:debug
```

### Python

```bash
# Run all tests
pytest

# Run specific test
pytest tests/unit/test_parser.py

# Run with coverage
pytest --cov=src

# Verbose output
pytest -v

# Debug mode
pytest --pdb
```

---

## Resources

- [Jest Documentation](https://jestjs.io/)
- [pytest Documentation](https://docs.pytest.org/)
- [VSCode Extension Testing](https://code.visualstudio.com/api/working-with-extensions/testing-extension)
- [Testing Best Practices](https://testingjavascript.com/)

---

## Questions?

- Open an issue with label `testing`
- Ask in GitHub Discussions
- Check existing test files for examples

---

**Last Updated**: 2026-02-28
**Version**: 1.0.0
