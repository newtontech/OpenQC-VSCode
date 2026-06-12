# Contributing to OpenQC-VSCode

Thank you for your interest in contributing to OpenQC-VSCode! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [How to Contribute](#how-to-contribute)
- [Pull Request Process](#pull-request-process)
- [Coding Standards](#coding-standards)
- [Testing Guidelines](#testing-guidelines)
- [Documentation](#documentation)

## Code of Conduct

### Our Pledge

We are committed to providing a welcoming and inspiring community for all. Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

### Our Standards

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive criticism
- Show empathy towards others

## Getting Started

### Prerequisites

- Node.js 18 or higher
- Python 3.8 or higher
- VSCode 1.85 or higher
- Git

### Fork and Clone

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/OpenQC-VSCode.git
   cd OpenQC-VSCode
   ```
3. Add upstream remote:
   ```bash
   git remote add upstream https://github.com/newtontech/OpenQC-VSCode.git
   ```

## Development Setup

### Install Dependencies

```bash
# Node.js dependencies
npm install

# Python dependencies
pip install -e .[dev]
```

### Build

```bash
# Compile TypeScript
npm run compile

# Watch mode
npm run watch
```

### Run Tests

```bash
# TypeScript tests
npm test

# Python tests
pytest

# With coverage
npm run test:coverage
pytest --cov=src
```

### Debug Extension

1. Open project in VSCode
2. Press `F5` to launch Extension Development Host
3. Test your changes in the new VSCode window

## How to Contribute

### Reporting Bugs

1. Check if the bug has already been reported in [Issues](https://github.com/newtontech/OpenQC-VSCode/issues)
2. If not, create a new issue using the **Bug Report** template
3. Include:
   - Clear description
   - Steps to reproduce
   - Expected vs actual behavior
   - Environment details
   - Screenshots (if applicable)

### Suggesting Features

1. Check existing [Issues](https://github.com/newtontech/OpenQC-VSCode/issues) for similar requests
2. Create a new issue using the **Feature Request** template
3. Describe:
   - Problem it solves
   - Proposed solution
   - Use cases
   - Alternatives considered

### Working on Issues

1. Look for issues labeled:
   - `good first issue` - Good for newcomers
   - `help wanted` - Community help needed
   - `bug` - Bug fixes
   - `feature` - New features

2. Comment on the issue to indicate you're working on it

3. Follow the [TDD Guidelines](docs/project/TDD-GUIDELINES.md) for development

## Pull Request Process

### 1. Create a Branch

```bash
# Update main branch
git checkout main
git pull upstream main

# Create feature branch
git checkout -b feature/issue-123-vasp-parser
```

### 2. Make Changes

- Follow [Coding Standards](#coding-standards)
- Write tests first (TDD)
- Update documentation
- Keep commits focused

### 3. Commit Changes

We follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat: add INCAR parser
fix: correct POSCAR coordinate parsing
docs: update installation instructions
test: add tests for Gaussian parser
refactor: simplify conversion logic
style: format code with prettier
chore: update dependencies
```

### 4. Push and Create PR

```bash
git push origin feature/issue-123-vasp-parser
```

Then create a Pull Request on GitHub.

### 5. PR Checklist

Before submitting, ensure:

- [ ] Tests pass (`npm test`, `pytest`)
- [ ] Coverage maintained (≥ 80%)
- [ ] Linter passes (`npm run lint`)
- [ ] Documentation updated
- [ ] CHANGELOG.md updated
- [ ] PR references issue (e.g., "Closes #123")
- [ ] PR template filled out

### 6. Code Review

- Respond to all comments
- Make requested changes
- Push new commits (don't force push)
- Request re-review when ready

### 7. Merge

A maintainer will merge your PR after:
- All CI checks pass
- At least one approval
- All comments resolved

## Coding Standards

### TypeScript

We use:
- **ESLint** for linting
- **Prettier** for formatting
- **TypeScript strict mode**

```bash
# Check linting
npm run lint

# Fix linting issues
npm run lint:fix

# Format code
npm run format
```

#### Code Style

```typescript
// Use interfaces for object shapes
interface ParseResult {
  success: boolean;
  data: unknown;
}

// Use async/await over promises
async function parseFile(path: string): Promise<ParseResult> {
  const content = await fs.readFile(path, 'utf-8');
  return parser.parse(content);
}

// Use descriptive names
const numberOfAtoms = atoms.length; // Good
const n = atoms.length; // Bad

// Document public APIs
/**
 * Parse VASP INCAR file
 * @param content - Raw file content
 * @returns Parsed parameters
 */
export function parseIncar(content: string): Record<string, unknown> {
  // ...
}
```

### Python

We follow PEP 8 with:
- **Black** for formatting
- **Flake8** for linting
- **MyPy** for type checking

```bash
# Format code
black src tests

# Check linting
flake8 src tests

# Type check
mypy src
```

#### Code Style

```python
from typing import Dict, List, Optional


def parse_incar(content: str) -> Dict[str, any]:
    """
    Parse VASP INCAR file.
    
    Args:
        content: Raw file content
        
    Returns:
        Dictionary of parameters
        
    Raises:
        ParseError: If content is invalid
    """
    # Use descriptive variable names
    parameters = {}
    
    # Document complex logic
    for line in content.split('\n'):
        # Skip comments and empty lines
        if line.startswith('#') or not line.strip():
            continue
            
        # Parse parameter
        key, value = parse_parameter(line)
        parameters[key] = value
    
    return parameters
```

## Testing Guidelines

We follow Test-Driven Development (TDD). See [TDD-GUIDELINES.md](docs/project/TDD-GUIDELINES.md) for details.

### Test Structure

```
tests/
├── unit/           # Fast, isolated tests
├── integration/    # Test component interactions
├── e2e/            # Test complete workflows
└── fixtures/       # Test data
```

### Writing Tests

```typescript
// tests/unit/parsers/incarParser.test.ts
import { IncarParser } from '@/parsers/vasp';

describe('IncarParser', () => {
  describe('parse()', () => {
    it('should parse basic parameters', () => {
      const input = 'ENCUT = 520\nISMEAR = 0';
      const result = IncarParser.parse(input);
      
      expect(result.success).toBe(true);
      expect(result.data).toEqual({
        ENCUT: 520,
        ISMEAR: 0
      });
    });

    it('should throw error on invalid input', () => {
      const input = 'INVALID = abc';
      expect(() => IncarParser.parse(input)).toThrow(ParseError);
    });
  });
});
```

### Coverage Requirements

- **Overall**: ≥ 80%
- **Critical paths** (parsers, converters): ≥ 90%
- **UI components**: ≥ 70%

## Documentation

### Code Documentation

- Document all public APIs
- Use JSDoc/TSDoc for TypeScript
- Use docstrings for Python
- Include examples

### User Documentation

- Update README.md for user-facing changes
- Update docs/ for detailed guides
- Add examples to docs/examples/

### Changelog

Update CHANGELOG.md:

```markdown
## [Unreleased]

### Added
- New INCAR parser with validation

### Changed
- Improved POSCAR parsing performance

### Fixed
- Bug in Gaussian coordinate parsing
```

## Project Structure

```
OpenQC-VSCode/
├── .github/              # GitHub configs
│   ├── ISSUE_TEMPLATE/   # Issue templates
│   └── workflows/        # CI/CD pipelines
├── docs/                 # Documentation
│   ├── architecture/     # Architecture docs
│   └── api/              # API reference
├── src/                  # TypeScript source
│   ├── parsers/          # File parsers
│   ├── converters/       # Format converters
│   ├── visualization/    # 3D rendering
│   └── ai/               # AI features
├── tests/                # Test suites
│   ├── unit/
│   ├── integration/
│   ├── e2e/
│   └── fixtures/
├── docs/project/PLAN.md          # Project roadmap
├── docs/project/TDD-GUIDELINES.md # Testing guide
├── docs/project/TASK-MANAGEMENT.md # Task management
└── README.md             # This file
```

## Getting Help

- 💬 [GitHub Discussions](https://github.com/newtontech/OpenQC-VSCode/discussions)
- 🐛 [GitHub Issues](https://github.com/newtontech/OpenQC-VSCode/issues)
- 📧 Email: support@newtontech.com
- 💬 Discord: [Join our community](https://discord.gg/openqc)

## Recognition

Contributors are recognized in:
- README.md contributors section
- Release notes
- GitHub contributors page

Thank you for contributing! 🎉


<!-- repo-governance-kit:contributing-v1 -->

## Agent-Ready Issue and PR Contract

Agent-ready work must have a GitHub issue with a goal, acceptance criteria, required tests, and out-of-scope notes. Branches should use `agent/issue-<number>-<slug>`, and PR bodies must include `Fixes #<issue-number>`.

Before requesting review, run `make format`, `make lint`, `make typecheck`, `make test`, and `make check`. Code changes should include tests or an explicit no-test justification in the PR body.
