# Contributing to CSA-3D

Thank you for your interest in contributing to CSA-3D! This document provides guidelines for contributing code, documentation, and bug reports.

## Code of Conduct

Please be respectful and constructive in all interactions with the community.

## Getting Started

### Setting Up Development Environment

1. Fork and clone the repository
2. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # on Windows: venv\Scripts\activate
   ```

3. Install in development mode:
   ```bash
   pip install -e ".[dev,docs]"
   ```

### Running Tests and Quality Checks

```bash
# Run tests
make test

# Run with coverage
make test-cov

# Check code style
make lint

# Auto-format code
make format

# Type checking
make type-check
```

## Development Workflow

1. Create a feature branch:
   ```bash
   git checkout -b feature/description
   ```

2. Make changes following code style guidelines (see below)

3. Add tests for new functionality

4. Ensure all tests pass:
   ```bash
   make lint
   make test
   make type-check
   ```

5. Commit with clear message:
   ```bash
   git commit -m "Fix: issue description

   - Detailed explanation if needed
   - Multiple points if necessary
   
   Fixes #123"
   ```

6. Push and create a pull request

## Code Style Guidelines

### Python Style

We follow PEP 8 with these tools:
- **black** for formatting (line length: 100)
- **isort** for import ordering
- **flake8** for linting

Run before committing:
```bash
make format  # Auto-fix
make lint    # Verify
```

### Type Hints

Add type hints to all public methods:
```python
def method(self, param: str, count: int = 0) -> bool:
    """Short description.
    
    Longer description.
    
    Args:
        param: Description of param
        count: Description of count. Defaults to 0.
    
    Returns:
        Description of return value
    
    Raises:
        ValueError: When validation fails
    """
```

### Docstrings

Use Google-style docstrings for all public classes and methods:
```python
class MyClass:
    """One-line summary.
    
    Longer description explaining the purpose, usage, and
    any important details.
    
    Attributes:
        attr1 (str): Description
        attr2 (int): Description
    
    Example:
        >>> obj = MyClass()
        >>> obj.method()
    """
    
    def method(self, param: str) -> str:
        """Short description.
        
        Args:
            param: Description
        
        Returns:
            Description of return value
        
        Raises:
            ValueError: When this happens
        """
```

### Naming Conventions

- `ClassName` - Classes (CamelCase)
- `function_name()` - Functions and methods (snake_case)
- `CONSTANT_NAME` - Constants (ALL_CAPS)
- `_private_method()` - Private methods (leading underscore)
- `__dunder__` - Special methods only

## Adding Tests

Tests are in the `tests/` directory using pytest.

### Test Structure

```python
import pytest
from csa3d import MyClass

@pytest.fixture
def my_object():
    """Fixture providing test data."""
    return MyClass(param="test")

def test_method_success(my_object):
    """Test successful case."""
    result = my_object.method()
    assert result == expected

def test_method_failure():
    """Test error handling."""
    obj = MyClass()
    with pytest.raises(ValueError):
        obj.invalid_method()

@pytest.mark.slow
def test_slow_operation():
    """Test that runs infrequently."""
    # ...
```

### Running Specific Tests

```bash
# Run single file
pytest tests/test_myfile.py

# Run single test
pytest tests/test_myfile.py::test_function

# Skip slow tests
pytest -m "not slow"

# Run with coverage
pytest --cov=csa3d tests/
```

## Documentation

### Updating README

- Keep examples simple and working
- Test examples before submitting
- Link to detailed docs for complex features

### Adding Documentation

- Docstrings in code (see examples above)
- Architecture guide: `docs/architecture.md`
- Usage examples: `examples/`
- API reference: Auto-generated from docstrings

### Building Docs

```bash
make docs
# Output in docs/_build/html/
```

## Reporting Issues

### Bug Reports

Include:
1. CSA-3D version
2. Python version
3. Operating system
4. Steps to reproduce
5. Expected vs actual behavior
6. Error message (full traceback)

Example:
```markdown
### Description
Template creation fails with multiple UnIProt references

### Steps to reproduce
1. Load entry 212
2. Call create_template()

### Error
```
NotImplementedError: Multiple UniProt references not supported
```

### Version info
- CSA-3D: 0.1.0
- Python: 3.10.4
- OS: macOS 12.4
```

### Feature Requests

Include:
1. Clear description of desired functionality
2. Motivation and use cases
3. Possible implementation approach
4. Related issues or discussions

Example:
```markdown
### Description
Add parallel processing for large batches

### Motivation
Building 100+ entries takes 30+ minutes sequentially

### Proposed approach
- Use ProcessPoolExecutor
- Each worker gets independent Mcsa instance
- Return collected results
```

## Pull Request Process

1. **Title**: Clear, concise description
   - ✓ `Fix: config path validation on non-existent dirs`
   - ✗ `Fixed bug`

2. **Description**: Explain changes
   - What problem does this solve?
   - How does it solve it?
   - Any breaking changes?
   - Related issues (#123)

3. **Checklist**:
   - [ ] Tests pass locally (`make test`)
   - [ ] Code style passes (`make lint`)
   - [ ] Type checking passes (`make type-check`)
   - [ ] Docstrings added/updated
   - [ ] README updated (if applicable)
   - [ ] No new warnings or errors

4. **Review**: Wait for maintainer feedback

## Performance Considerations

When modifying core classes:
1. Profile before/after: `pytest --durations=10 tests/`
2. Check memory usage for large datasets
3. Document performance impact in docstrings
4. Consider lazy loading for expensive operations

## Release Process

Only maintainers:
1. Update version in `pyproject.toml` and `csa3d_version.py`
2. Update `CHANGELOG.md`
3. Create release branch: `release/v0.2.0`
4. Tag release: `git tag v0.2.0`
5. Push and create release on GitHub
6. Build and upload to PyPI: `python -m build && python -m twine upload dist/*`

## Getting Help

- **Questions**: Open a discussion or issue
- **Bug help**: Comment on the issue
- **Design discussion**: Open an issue with `[DISCUSSION]` prefix
- **Chat**: TBD (Slack, Discord, etc.)

## Recognition

Contributors will be:
1. Listed in CONTRIBUTORS.md
2. Acknowledged in release notes
3. Thanked in README

## License

By contributing, you agree that your code will be licensed under Apache 2.0 (same as CSA-3D).

---

Thank you for making CSA-3D better! 🎉
