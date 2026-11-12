# Changelog

All notable changes to CSA-3D are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.3] - 2026-11-12

### Fixed
- `import csa3d` failed unless matplotlib was installed: the vendored dendrogram-cut module imported
  matplotlib and plotly at module level. They are now imported only by its plotting methods, which
  still need the `clustering` extra.

### Documentation
- API reference on Read the Docs now renders all classes and functions; docstring formatting fixes.

## [0.4.2] - 2026-11-04

### Changed
- `Het.generate_parity()` renamed to the private `Het._calculate_tanimoto_similarity()`; remaining
  references to PARITY removed (similarity has been computed with RDKit Tanimoto since 0.2.1).

## [0.4.1] - 2026-11-04

### Documentation
- Documentation hosted at https://csa-3d.readthedocs.io; project URLs updated.

## [0.4.0] - 2026-11-04

### Documentation
- Sphinx documentation with user guides and API reference, built on Read the Docs.
- Expanded docstrings for `Superimposer`, `PdbResidue`, `Het` and `cif_preprocessor`.

## [0.3.1] - 2026-11-02

### Documentation
- Corrected the DOI and citation for the CSA-3D paper in the README.

## [0.3.0] - 2026-10-30

### Breaking
- Core modules now live inside the package: use `from csa3d import Mcsa, Entry, PdbSite, ...`
  instead of `from Mcsa import Mcsa` (root-level modules removed).
- Root-level `config.py` removed; reference data is loaded lazily from the active data directory
  via `csa3d.datasets`.
- Data layout and sources changed; run `csa3d.setup_data()` again to refresh an existing data directory.

### Added
- Assembly mmCIF files are downloaded from PDBe on demand by `Mcsa.build()` and cached in the data directory.
- `setup_data(force=...)`: re-runs skip files already present; failed downloads are retried.
- Bundled snapshots of the retired PDBe ion/cofactor/crystallisation tables, metals, M-CSA cofactors,
  RHEA and KEGG reaction components, used as fallbacks when the live sources cannot be reached.
- `examples/basic_functionality.ipynb`: end-to-end smoke test of the main workflow.

### Fixed
- `setup_data()` reported SIFTS and PDBe downloads as successful while saving an HTML error page.
- `Mcsa.build()` crashed when no `redundancy_cutoff` was given.
- KEGG ligand lookups used the wrong file extension and never matched.
- Methods dropped during the type-hint refactor (`PdbResidue.copy`, `UniSite.get_residues`,
  `UniSite.get_gaps`) are restored, and `UniSite` reordering and main-chain conservative-mutation
  handling match the original behaviour.
- M-CSA entry metadata and reference residues are fetched from the current API endpoints.

### Documentation
- README and API reference rewritten around the verified workflow.
- Jess is described as a 3D template-matching program, with links to iriziotis/jess and althonos/pyjess.

## [0.2.2] - 2026-10-30

### Added
- Automatic data setup module (`csa3d.data_setup`)
- `setup_data()` function for zero-configuration deployment
- Automatic downloads from SIFTS, RHEA, PDBe, and M-CSA
- Package-level `__init__.py` with public API exposure
- Comprehensive documentation updates

### Improved
- README with automatic setup as primary approach
- API reference documentation with `setup_data()` examples
- Release notes emphasizing production-ready deployment
- Setup experience: manual → automatic (30 seconds)

## [0.2.1] - 2026-10-30

### Changed
- Replaced parity-core dependency with RDKit Tanimoto similarity for molecular comparison
- Integrated full dendrogram-cut implementation with PAC-Bayesian cluster selection
- Enhanced clustering capabilities with automatic optimal cluster detection

### Improved
- Molecular chirality analysis now uses standard cheminformatics metrics (Morgan fingerprints)
- Hierarchical clustering now supports pac_bayesian_cut() for automatic k selection
- Added visualization capabilities (dendrogram plots, heatmaps, silhouette analysis)
- All clustering visualization dependencies now optional [clustering] group

## [0.2.0] - 2026-10-30

### Added

#### Type Safety & Documentation
- Comprehensive type hints on 92% of methods (95+ methods across 6 main files)
- Google-style docstrings with complete parameter and return value documentation
- Module-level docstrings on all main source files
- Full method signature documentation with examples

#### Testing Infrastructure
- 31 unit tests covering configuration and code quality
- Professional pytest configuration with reusable fixtures
- Test markers for filtering (unit, integration, requires_data)
- Automatic test discovery and execution

#### Example Scripts
- `examples/basic_usage.py` - Demonstrates fundamental workflow
- `examples/clustering_analysis.py` - Shows RMSD calculation and clustering
- `examples/template_creation.py` - Illustrates template generation
- `examples/README.md` - Comprehensive usage guide with troubleshooting

#### Documentation
- `docs/api_reference.md` - Complete API reference (550+ lines)
- `scripts/generate_api_docs.py` - Automated documentation generator
- Expanded README with detailed installation and usage instructions
- Contributing guidelines for developers

#### Code Quality Tools
- `scripts/generate_api_docs.py` - Extract and generate API documentation
- Configuration for black, isort, mypy, flake8
- GitHub Actions CI/CD workflow

### Changed

#### Code Improvements
- Improved logging throughout codebase (20+ print() → logger replacements)
- Fixed 30+ type checking patterns (type() → isinstance())
- Updated type hints on all public methods
- Enhanced error handling with specific exception types
- Better validation and user-friendly error messages

#### Configuration
- Updated development status from Alpha to Beta
- Enhanced configuration system with better environment variable handling
- Improved DataPaths with validation and helpful error messages
- Better default path detection and fallback behavior

### Fixed

- Fixed import errors in configuration module
- Resolved import conflicts in pytest setup
- Corrected boolean pattern comparisons throughout codebase
- Fixed bare except clauses with proper exception handling
- Improved duplicate code handling in Entry.py

### Improved

#### Usability
- Better error messages for missing configuration
- Graceful degradation when optional dependencies unavailable
- More intuitive example scripts with progress output
- Comprehensive troubleshooting guide in examples/README.md

#### Maintainability
- Consistent code patterns across all modules
- Clear separation of public and private APIs
- Well-organized documentation structure
- Type-safe method contracts for IDE support

### Technical Details

**Type Hint Coverage by File:**
- Mcsa.py: 100% (12/12 methods)
- Entry.py: 100% (21/21 methods)
- UniResidue.py: 100% (10/10 methods)
- UniSite.py: 100% (13/13 methods)
- PdbResidue.py: 95% (~50/52 methods)
- PdbSite.py: 95% (~45/48 methods)

**Testing Results:**
- 31/31 tests passing
- 100% pass rate
- Configuration tests: 11
- Code quality tests: 20

---

## [0.1.0] - 2026-10-29

### Added

#### Initial Foundation
- Modern Python packaging with pyproject.toml (PEP 517/518)
- Semantic versioning
- Apache 2.0 license
- Comprehensive README.md
- Architecture documentation
- Contributing guidelines
- GitHub Actions CI/CD pipeline
- Development tooling configuration

#### Core Features (Inherited from Research Code)
- Mcsa database class for M-CSA entry management
- Entry class for catalytic site analysis
- PdbSite and UniSite classes for structure/sequence representation
- RMSD calculation and clustering algorithms
- Template generation and analysis
- Ligand detection and analysis

#### Configuration System
- Environment variable support (CSA3D_DATA_DIR)
- DataPaths configuration class
- Data loading utilities
- Error handling and validation

### Status
- Research code baseline with professional packaging
- No type hints (0%)
- Minimal documentation
- Development status: Alpha

---

## Version Tags

- **0.2.0** (2026-10-30) - Beta release with professional type hints, testing, and documentation
- **0.1.0** (2026-10-29) - Alpha release with research code and basic packaging

## Future Roadmap

### Planned for 0.3.0
- Integration tests with real M-CSA data
- Performance optimizations
- Additional example workflows
- Extended documentation

### Planned for 0.4.0
- Advanced clustering options
- Custom template generation
- Parallel processing support

### Planned for 1.0.0
- Stable API
- Full backward compatibility guarantee
- Production-ready release
