# CSA-3D Refactoring Plan

## Overview
This document outlines a comprehensive refactoring to improve usability, documentation, maintainability, and automatic data setup for the CSA-3D package.

---

## Phase 1: Foundation (Setup & Configuration)

### 1.1 Package Structure & Setup Files
**Files to create:**
- `pyproject.toml` - Modern Python packaging (replaces setup.py)
- `requirements.txt` - Core dependencies
- `requirements-dev.txt` - Development dependencies
- `.python-version` - Python version specification
- `MANIFEST.in` - Include non-Python files
- `LICENSE` - License file (recommend Apache 2.0 or MIT for public use)

**Key decisions:**
- Use Python 3.8+ as minimum (reasonable for scientific packages)
- Pin major versions, allow minor updates
- Separate core, dev, and optional dependencies

### 1.2 Configuration Management System
**Replace `config.py` with new `csa3d/config/`:**

```
csa3d/
├── config/
│   ├── __init__.py
│   ├── settings.py      # Configuration class with env var support
│   ├── data_paths.py    # Path management
│   └── loader.py        # Data loading with caching and error handling
```

**Implementation:**
- Load paths from environment variables first, then fallback to defaults
- Support config file (YAML/TOML) for user overrides
- Validate required files on import, give helpful errors
- Lazy-load data (only when accessed)
- Add `setup_data()` function to download missing files

### 1.3 Data Download Script
**New file: `scripts/download_data.py`**
- Script to fetch required datasets from EMBL-EBI, PDBe, KEGG, RHEA, SIFTS
- Create directory structure automatically
- Resume capability for large files
- Verification (checksums if available)
- Progress reporting

---

## Phase 2: Documentation

### 2.1 README Complete Rewrite
**Sections:**
- What is CSA-3D (motivation, use cases)
- Installation (pip, from source)
- Quick start (minimal working example)
- Setup (data requirements, download)
- Main usage patterns
- Examples with output
- Contributing guidelines
- Citation info (from your paper)

### 2.2 API Documentation
**Add docstrings to all classes/methods:**

Template:
```python
def method(self, param1: Type1, param2: Type2 = default) -> ReturnType:
    """Short description.
    
    Longer description explaining what it does and why.
    
    Args:
        param1: Description of param1
        param2: Description of param2. Defaults to default.
    
    Returns:
        Description of return value
    
    Raises:
        ExceptionType: When this exception occurs
    
    Examples:
        >>> obj.method(1, 2)
        expected_output
    """
```

**Priority order:**
1. `Mcsa` class - main entry point
2. `Entry` class - core data structure
3. `PdbSite`, `UniSite` classes
4. Other classes and functions

### 2.3 Architecture Documentation
**New file: `docs/architecture.md`**
- Data flow diagram (ASCII art)
- Class relationships
- Key concepts (reference structures, homologues, etc.)
- Design decisions and trade-offs

### 2.4 Usage Examples
**New directory: `examples/`**
- `basic_usage.py` - Load entry, inspect sites
- `template_creation.py` - Create and write templates
- `clustering_analysis.py` - RMSD clustering
- `batch_processing.py` - Process multiple entries
- `jupyter_notebook.ipynb` - Interactive walkthrough

---

## Phase 3: Code Quality

### 3.1 Type Hints Throughout
**Priority:**
1. Public APIs (method signatures)
2. Class attributes (in `__init__`)
3. Internal methods
4. Complex functions with many parameters

**Tools:** Use `pyright` or `mypy` for validation

### 3.2 Code Quality Fixes
**Anti-patterns to fix:**
- `type(x) == Class` → `isinstance(x, Class)`
- Add `__all__` to `__init__.py`
- Replace string type checks with actual type hints
- Consistent error handling patterns
- Use logging instead of print statements (optional: for warnings)

**Code style:**
- Run `black` for consistent formatting
- Run `isort` for import ordering
- Add `flake8` config to catch issues

### 3.3 Error Handling
**Improvements:**
- Validate input parameters early
- Provide specific error messages (not silent failures)
- Add logging for debugging
- Graceful degradation where appropriate
- Custom exceptions for library-specific errors

---

## Phase 4: Testing

### 4.1 Unit Tests Structure
**Create `tests/` directory:**
```
tests/
├── conftest.py              # Pytest fixtures
├── test_config.py           # Config loading/validation
├── test_mcsa.py             # Mcsa class
├── test_entry.py            # Entry class
├── test_pdbsite.py          # PdbSite class
├── test_unisite.py          # UniSite class
└── fixtures/                # Test data
    └── minimal_data/        # Small test datasets
```

**Approach:**
- Start with unit tests for config and data loading
- Test core Entry operations
- Use fixtures for test data (don't require full DB)
- Aim for >70% coverage initially, >85% eventually
- Mark slow tests with `@pytest.mark.slow`

### 4.2 CI/CD Configuration
**New files:**
- `.github/workflows/tests.yml` - Run tests on push/PR
- `.github/workflows/linting.yml` - Code quality checks
- `tox.ini` - Test matrix (Python 3.8, 3.9, 3.10, 3.11)

---

## Phase 5: Usability Improvements

### 5.1 Public API Cleanup
**`__init__.py` should export:**
```python
from .Mcsa import Mcsa
from .Entry import Entry
from .PdbSite import PdbSite
from .UniSite import UniSite

__all__ = ['Mcsa', 'Entry', 'PdbSite', 'UniSite']
```

**Consider deprecating/hiding:**
- Direct use of `PdbResidue`, `UniResidue` (internal implementation detail)
- Private methods (prefix with `_`)

### 5.2 Configuration API
**Add convenience functions:**
```python
csa3d.config.setup_data(data_dir)  # Download/configure
csa3d.config.validate()            # Check if setup complete
csa3d.config.get_data_dir()        # Get current config
```

### 5.3 Better Error Messages
**Examples:**
```
Current:  "Could not find mmCIF file for 1abc"
Better:   "Could not find mmCIF file for 1abc. Run `csa3d.config.setup_data()` 
          to download required data files."

Current:  "Attempted to add non-Site object in entry"
Better:   "Cannot add object of type {type(site).__name__} to entry. 
          Expected PdbSite or UniSite."
```

---

## Phase 6: Additional Improvements

### 6.1 Caching & Performance
- Implement caching for parsed mmCIF files
- Cache clustering results
- Document performance characteristics

### 6.2 Optional Dependencies
- Make visualization (if any) optional
- Document all external requirements
- Handle gracefully if optional packages missing

### 6.3 Version Management
- Semantic versioning
- Version file in `csa3d/__version__.py`
- Update in multiple places with CI check

---

## Implementation Order

### Sprint 1: Foundation
- [ ] Create pyproject.toml
- [ ] Create requirements.txt files
- [ ] Refactor config.py → config/ package
- [ ] Add environment variable support
- [ ] Add basic validation and error messages

### Sprint 2: Documentation
- [ ] Complete README
- [ ] Add docstrings to Mcsa, Entry, PdbSite, UniSite
- [ ] Create architecture.md
- [ ] Create basic_usage.py example

### Sprint 3: Code Quality
- [ ] Add type hints to public APIs
- [ ] Fix anti-patterns (type checks, etc.)
- [ ] Add code style checks (black, isort, flake8)
- [ ] Improve error messages

### Sprint 4: Testing & Data
- [ ] Create conftest.py with fixtures
- [ ] Create minimal test dataset
- [ ] Write unit tests for config and core classes
- [ ] Create data download script

### Sprint 5: Polish
- [ ] Complete remaining docstrings
- [ ] Add more examples
- [ ] Create CI/CD workflows
- [ ] Documentation website (optional)

---

## File Changes Summary

### New Files
- `pyproject.toml`
- `requirements.txt`, `requirements-dev.txt`
- `csa3d/config/` (package)
- `csa3d/__version__.py`
- `docs/architecture.md`
- `examples/` (multiple files)
- `tests/` (multiple files)
- `scripts/download_data.py`
- `.github/workflows/` (CI/CD)
- `LICENSE`

### Modified Files
- `csa3d/__init__.py` - Better exports
- `csa3d/config.py` → Migrate to `csa3d/config/`
- All class files - Add type hints and docstrings
- `README.md` - Complete rewrite

### Deleted/Refactored
- Move `config.py` → `config/settings.py`
- Move data loading logic → `config/loader.py`

---

## Success Criteria

1. **Installation**: `pip install csa3d` works on fresh machine
2. **Setup**: Running `csa3d.config.setup_data()` downloads required files
3. **Documentation**: Every public class/method has examples
4. **Testing**: >70% code coverage with passing tests
5. **Usability**: New user can load entry in <5 lines of code
6. **Quality**: Passes black, isort, flake8, mypy checks
7. **CI/CD**: Tests pass on Python 3.8+

