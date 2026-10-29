# CSA-3D Architecture

## System Overview

CSA-3D is designed to work with catalytic site data from three primary sources:

1. **M-CSA (Mechanism and Catalytic Site Atlas)** - Reference catalytic sites
2. **PDB** - 3D protein structures from crystallography
3. **UniProt** - Protein sequences and annotations

The package bridges these data sources to enable structural analysis and comparison of enzyme active sites.

## Data Flow

```
┌─────────────────────────────────────────────────────────────┐
│                    External Data Sources                     │
├──────────────────────┬──────────────────┬──────────────────┤
│    M-CSA JSON        │   PDB mmCIF      │  UniProt Seq     │
│  (Catalytic          │  (Structures)    │  (Annotations)   │
│   Residues)          │                  │                  │
└──────────────────────┴──────────────────┴──────────────────┘
           │                   │                   │
           └───────────────────┼───────────────────┘
                               │
                               ▼
                ┌──────────────────────────────┐
                │   csa3d.config (Settings)    │
                │  • Path management           │
                │  • Data validation           │
                │  • Lazy loading              │
                └──────────────────────────────┘
                               │
                ┌──────────────┴──────────────┐
                │                             │
                ▼                             ▼
         ┌────────────────────┐      ┌──────────────────┐
         │ PdbResidue (Atom)  │      │ UniResidue       │
         │ Het (Ligand)       │      │ (Sequence Pos.)  │
         └────────────────────┘      └──────────────────┘
                │                             │
                ▼                             ▼
         ┌────────────────────┐      ┌──────────────────┐
         │  PdbSite           │      │  UniSite         │
         │ • 3D Coordinates   │      │ • Seq Alignment  │
         │ • Geometry         │      │ • Orthologs      │
         │ • Ligands          │      │ • Annotations    │
         └────────────────────┘      └──────────────────┘
                │                             │
                └──────────────┬──────────────┘
                               │
                               ▼
                      ┌────────────────────┐
                      │    Entry           │
                      │ • Sites (PDB+Uni)  │
                      │ • Clustering       │
                      │ • Templates        │
                      └────────────────────┘
                               │
                               ▼
                      ┌────────────────────┐
                      │    Mcsa (Db)       │
                      │ • All Entries      │
                      │ • Batch Loading    │
                      │ • Iteration        │
                      └────────────────────┘
```

## Core Classes

### Configuration Layer

**`csa3d.config.settings`** - Path and configuration management
- `DataPaths` - Manages file paths, environment variables
- `ConfigError` - Configuration-related errors

**`csa3d.config.loader`** - Data loading utilities
- Load functions for each data type (JSON, CSV)
- Lazy loading with error handling
- Caching support

### Data Layer

**`csa3d.PdbResidue`** - Represents a catalytic residue in a PDB structure
- Attributes: `pdb_id`, `chain`, `resnum`, `resname`
- Methods: `add_structure()`, `get_distance()`, `is_standard`
- Relationships: Contains Bio.PDB `Residue` structure

**`csa3d.PdbResidue.Het`** - Represents a ligand/heteroatom in PDB
- Attributes: `resname`, `type`, `similarity`
- Methods: `get_distance()` to residues

**`csa3d.UniResidue`** - Represents a catalytic residue position in UniProt
- Attributes: `uniprot_id`, `position`, `resname`
- Methods: `get_alignment()`, matching to PdbResidue

### Site Layer

**`csa3d.PdbSite`** - A complete catalytic site from a PDB structure
- Attributes: `residues`, `ligands`, `pdb_id`, `parent_structure`
- Methods:
  - `build_reference()` - Create from reference residues
  - `build_all()` - Create from homologous residues
  - `fit()` - Structural superposition with another site
  - `find_ligands()` - Identify proximal ligands
  - `get_functional_site()` - Extract functional atoms only

**`csa3d.UniSite`** - A complete catalytic site from UniProt sequence
- Attributes: `uniprot_id`, `residues`
- Methods:
  - `build()` - Create from residue list
  - Sequence-level comparisons

### Entry Layer

**`csa3d.Entry`** - Represents one M-CSA entry
- Attributes:
  - `mcsa_id` - M-CSA identifier
  - `pdbsites` - List of all PdbSite objects
  - `unisites` - List of all UniSite objects
  - `info` - Metadata (EC number, enzyme name, etc.)

- Key Methods:
  - `rmsd_matrix()` - Calculate pairwise RMSD between sites
  - `clustering()` - Hierarchical clustering by RMSD
  - `clustering_bayesian()` - Automatic cluster detection
  - `create_template()` - Generate active site template
  - `break_template()` - Split template into residue groups
  - `write_template()` - Export in PDB format

### Database Layer

**`csa3d.Mcsa`** - Top-level database class
- Attributes:
  - `entries` - Dictionary of Entry objects
  - `catalytic_residue_info` - M-CSA JSON data
  - `mcsa_entry_info` - Entry metadata

- Key Methods:
  - `build()` - Load and build Entry objects from M-CSA
  - `add()` - Add Entry to database
  - `update_pdbs()` - Download PDB structures
  - `__iter__()` - Iterate over all entries

## Key Data Structures

### Residue Identity

Residues are identified by:
```python
full_id = (pdb_id, chain, resnum, insertion_code, resname)
```

Example: `('1abc', 'A', 42, ' ', 'HIS')`

### Site Matching

Sites are matched by:
1. **Residue-by-residue alignment** - Based on SIFTS mapping
2. **Structural similarity** - RMSD threshold
3. **Sequence conservation** - Conservation vs. mutation status

### Ligand Classification

Ligands are classified by:
- **Type**: Substrate, cofactor, ion, artefact
- **Similarity**: Chemical similarity to cognate ligand (0-1)
- **Distance**: Minimum distance to any catalytic residue

## Workflow Examples

### Load and Build Entry

```
Mcsa.build(entry_ids=1)
  └─> Entry(mcsa_id=1) created
      └─> _build_pdb_residues()
          └─> PdbResidue.from_json() x N
      └─> _build_uniprot_residues()
          └─> UniResidue.from_json() x N
      └─> _build_pdb_sites()
          └─> PdbSite.build_reference() → reference site
          └─> PdbSite.build_all() → homologous sites
      └─> _build_uniprot_sites()
          └─> UniSite.build() x N
```

### Clustering and Template Creation

```
Entry.rmsd_matrix(sites)
  └─> Pairwise site.fit() → RMSD matrix

Entry.clustering(rmsd_matrix)
  └─> Hierarchical linkage
  └─> Cut dendrogram → clusters

Entry.create_template()
  └─> Select reference site
  └─> Superimpose all conserved sites
  └─> Calculate average coordinates
  └─> Find template representative
  └─> Write PDB with constraints
```

## Performance Characteristics

| Operation | Time | Memory | Notes |
|-----------|------|--------|-------|
| Build entry | 10-60s | 100-500 MB | Depends on # sites and PDB size |
| RMSD matrix (20 sites) | 30-120s | 100 MB | Quadratic in site count |
| Clustering (20 sites) | <1s | 50 MB | Fast after RMSD computed |
| Template creation | 5-30s | 100-300 MB | Fitting all sites |
| Template writing | <1s | 10 MB | Output file only |

## Design Decisions

### Why Type Hints Were Missing
Original development focused on functionality over code quality. Type hints are being added incrementally.

### Why Large Files
Classes like `PdbSite` (1046 lines) handle multiple concerns:
- Structure parsing (Bio.PDB integration)
- Geometry calculations
- Ligand detection
- File output

Future: Split into separate modules (parsing, geometry, ligand detection, output).

### Why Lazy Loading
Data files are large (GB). Lazy loading allows importing package without full data load.
Cost: More complex initialization flow.

### Clustering Algorithm
Uses hierarchical linkage with automatic height detection (70% of max distance).
Alternative: Bayesian cut available for probabilistic selection.

## Extension Points

### Adding New Entry Types
1. Subclass `Entry` or create new site type
2. Implement `build_from_data()` method
3. Add to `Mcsa._build_sites()` dispatcher

### Adding New Ligand Types
1. Define in `residue_definitions.py`
2. Update ligand detection in `PdbSite.find_ligands()`

### Custom Clustering
1. Subclass `Entry` or use as-is
2. Override `clustering()` method
3. Return `{cluster_id: [site_ids]}`

## External Dependencies

**Core Scientific Libraries:**
- `biopython` - PDB/mmCIF parsing, structure handling
- `numpy` - Numerical operations, arrays
- `scipy` - Clustering, distance metrics
- `pandas` - Tabular data

**Specialized Tools:**
- `rdkit` - Molecule fingerprints and similarity
- `scikit-learn` - K-means clustering
- `rmsd` - Optimal structural superposition
- `dendrogram-cut` - Bayesian cluster detection
- `parity_core` - (Used in original, purpose unclear)

## Testing Strategy

**Unit Tests** (high coverage target):
- Config loading and validation
- Residue/site creation from JSON
- RMSD calculations
- Template generation

**Integration Tests** (slower, require data):
- Full entry build
- Multi-site clustering
- File writing and parsing

**Fixtures:**
- Minimal test data (< 1 MB)
- Cached for fast test runs
- Not included in git (generated per test)

---

See README.md for usage examples and API documentation.
