# Contributing to OpenQC

Thank you for your interest in contributing to OpenQC! This document provides guidelines and instructions for contributing.

## 🌟 Ways to Contribute

- **Bug Reports**: Submit issues for bugs you encounter
- **Feature Requests**: Suggest new features or improvements
- **Code Contributions**: Submit pull requests
- **Documentation**: Improve or add documentation
- **Examples**: Add example input files or workflows

## 🛠️ Development Setup

### Prerequisites

- Python 3.9+
- Node.js 18+
- VSCode

### Setup Steps

1. **Fork and Clone**
   ```bash
   git clone https://github.com/YOUR_USERNAME/OpenQC.git
   cd OpenQC
   ```

2. **Python Development**
   ```bash
   cd core
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -e ".[dev]"
   ```

3. **VSCode Extension Development**
   ```bash
   cd ../vscode-extension
   npm install
   npm run watch
   ```

4. **Run Tests**
   ```bash
   # Python tests
   cd core
   pytest

   # Extension tests
   cd ../vscode-extension
   npm test
   ```

## 📝 Coding Standards

### Python

- Follow PEP 8 style guide
- Use type hints
- Write docstrings for all functions
- Maximum line length: 100 characters

```python
def parse_file(filepath: Path) -> ParsedInput:
    """Parse a quantum chemistry input file.
    
    Args:
        filepath: Path to the input file
        
    Returns:
        ParsedInput object with structure and parameters
    """
    ...
```

### TypeScript

- Use TypeScript strict mode
- Follow VSCode extension guidelines
- Document all public functions

## 🧪 Testing

### Running Tests

```bash
# All tests
pytest

# With coverage
pytest --cov=openqc --cov-report=html

# Specific test file
pytest tests/test_gaussian_parser.py
```

### Writing Tests

- Place tests in `tests/` directory
- Name test files `test_*.py`
- Use descriptive test names

```python
def test_gaussian_parser_parses_basic_input():
    """Test that GaussianParser correctly parses a basic input file."""
    parser = GaussianParser()
    result = parser.parse(SAMPLE_INPUT)
    assert result.structure is not None
    assert result.structure.num_atoms == 3
```

## 📚 Adding New Parsers

1. Create a new parser in `core/openqc/parsers/`
2. Inherit from `BaseParser`
3. Implement the `parse` method
4. Add tests
5. Update documentation

```python
# core/openqc/parsers/myformat.py
from openqc.parsers.base import BaseParser, ParsedInput

class MyFormatParser(BaseParser):
    def parse(self, content: str) -> ParsedInput:
        # Parse the content
        ...
        return ParsedInput(...)
```

## 🔧 Adding MCP Tools

1. Add tool definition to `list_tools()` in `ai-protocols/mcp-server/server.py`
2. Implement handler in `call_tool()`
3. Add tests
4. Update documentation

## 📖 Documentation

- Update README.md for user-facing changes
- Add docstrings for new functions
- Update API documentation for new features
- Add examples for complex features

## 🐛 Bug Reports

When submitting a bug report, please include:

1. OpenQC version
2. Python/Node.js version
3. Operating system
4. Steps to reproduce
5. Expected behavior
6. Actual behavior
7. Input file (if applicable, sanitized)

## 🎯 Pull Request Process

1. Create a feature branch
   ```bash
   git checkout -b feature/my-feature
   ```

2. Make your changes
   - Follow coding standards
   - Add tests
   - Update documentation

3. Run tests and linting
   ```bash
   pytest
   npm run lint
   ```

4. Commit with descriptive message
   ```bash
   git commit -m "feat: add support for CP2K input files"
   ```

5. Push and create PR
   ```bash
   git push origin feature/my-feature
   ```

6. Wait for review and address feedback

## 📋 Commit Message Format

Follow conventional commits:

- `feat:` New feature
- `fix:` Bug fix
- `docs:` Documentation changes
- `style:` Code style changes (formatting)
- `refactor:` Code refactoring
- `test:` Adding tests
- `chore:` Maintenance tasks

Examples:
```
feat: add CP2K input file parser
fix: correct lattice vector parsing in VASP parser
docs: update installation instructions
test: add tests for ORCA parser
```

## 🤝 Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Help others learn and grow

## 📞 Getting Help

- Open an issue for questions
- Join discussions in existing issues
- Check documentation first

Thank you for contributing! 🎉
