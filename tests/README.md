# Test Suite - OpenQC-VSCode

This directory contains all tests for the OpenQC-VSCode extension.

## Structure

```
tests/
├── unit/           # Fast, isolated unit tests
├── integration/    # Integration tests with real components
├── e2e/            # End-to-end tests (VSCode UI)
└── fixtures/       # Test data and sample files
    ├── vasp/       # VASP input file samples
    ├── gaussian/   # Gaussian input file samples
    └── orca/       # ORCA input file samples
```

## Running Tests

### TypeScript Tests

```bash
# Run all tests
npm test

# Run unit tests only
npm run test:unit

# Run integration tests
npm run test:integration

# Run e2e tests
npm run test:e2e

# Run with coverage
npm run test:coverage

# Watch mode
npm run test:watch
```

### Python Tests

```bash
# Run all tests
pytest

# Run unit tests
pytest tests/unit

# Run with coverage
pytest --cov=src

# Verbose output
pytest -v

# Specific test file
pytest tests/unit/test_parser.py
```

## Test Categories

### Unit Tests (`unit/`)
- Fast execution (< 100ms each)
- No external dependencies
- Test individual functions/classes in isolation
- Mock all external calls

### Integration Tests (`integration/`)
- Test interactions between components
- May use real file system or VSCode APIs
- Slower but more realistic
- Run after unit tests pass

### E2E Tests (`e2e/`)
- Test complete user workflows
- Use VSCode testing framework
- Run in CI only (not locally)
- Slowest but most comprehensive

## Writing Tests

Follow the [TDD Guidelines](../TDD-GUIDELINES.md) for best practices.

### Example Unit Test

```typescript
// tests/unit/parsers/vasp/incarParser.test.ts
import { IncarParser } from '@/parsers/vasp';

describe('IncarParser', () => {
  it('should parse basic INCAR parameters', () => {
    const input = 'ENCUT = 520\nISMEAR = 0';
    const result = IncarParser.parse(input);
    expect(result).toEqual({ ENCUT: 520, ISMEAR: 0 });
  });
});
```

### Example Integration Test

```python
# tests/integration/test_conversion.py
import pytest
from openqc.converters import VASPConverter

def test_convert_vasp_to_gaussian(tmp_path):
    """Test full conversion workflow."""
    poscar = tmp_path / "POSCAR"
    poscar.write_text(SAMPLE_POSCAR)
    
    converter = VASPConverter()
    result = converter.to_gaussian(str(poscar))
    
    assert result.success
    assert (tmp_path / "input.com").exists()
```

## Test Fixtures

Sample files are stored in `tests/fixtures/`:

- **vasp/** - POSCAR, INCAR, KPOINTS, POTCAR samples
- **gaussian/** - .com, .gjf input files
- **orca/** - .inp input files

Use these fixtures in tests:

```typescript
import { readFileSync } from 'fs';
import { join } from 'path';

const fixture = readFileSync(
  join(__dirname, '../fixtures/vasp/POSCAR-Si'),
  'utf-8'
);
```

## Coverage Requirements

- **Overall**: ≥ 80%
- **Critical paths** (parsers, converters): ≥ 90%
- **UI components**: ≥ 70%

Check coverage:
```bash
npm run test:coverage
# Open coverage/lcov-report/index.html
```

## CI/CD

Tests run automatically on:
- Every push to main/master/develop
- Every pull request
- Before releases

See `.github/workflows/test.yml` for details.

## Debugging Tests

### TypeScript

```bash
# Debug in VSCode
npm run test:debug

# Then use VSCode debugger (F5)
```

### Python

```bash
# Drop into debugger on failure
pytest --pdb

# Debug specific test
pytest tests/unit/test_parser.py::test_function -v
```

## Best Practices

1. ✅ Write tests first (TDD)
2. ✅ Keep tests focused and independent
3. ✅ Use descriptive test names
4. ✅ Test edge cases and errors
5. ✅ Keep fixtures minimal but realistic
6. ✅ Clean up after tests
7. ✅ Run tests before committing

## Resources

- [Jest Documentation](https://jestjs.io/)
- [pytest Documentation](https://docs.pytest.org/)
- [VSCode Extension Testing](https://code.visualstudio.com/api/working-with-extensions/testing-extension)
