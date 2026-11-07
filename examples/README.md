# CSA-3D Example Scripts

This directory contains example scripts demonstrating various CSA-3D workflows.

## Requirements

Before running examples, ensure CSA-3D data is available:

```bash
export CSA3D_DATA_DIR=/path/to/csa3d_data
```

Or place data in the default location: `~/.csa3d/data`

## Examples

### 1. Basic Usage (`basic_usage.py`)

Demonstrates the fundamental workflow:
- Initializing the M-CSA database
- Building a specific entry
- Accessing catalytic sites
- Calculating RMSD between sites

**Usage:**
```bash
python examples/basic_usage.py
```

**Output:**
- Database initialization status
- Entry information (EC number, number of sites)
- Site details (PDB ID, number of residues)
- RMSD matrix statistics

### 2. Clustering Analysis (`clustering_analysis.py`)

Demonstrates structural analysis of catalytic sites:
- Loading an entry
- Filtering conserved sites
- Calculating RMSD matrix
- Performing hierarchical clustering
- Analyzing cluster composition

**Usage:**
```bash
python examples/clustering_analysis.py
```

**Output:**
- RMSD statistics (min, max, mean distances)
- Number of clusters identified
- Cluster compositions

### 3. Template Creation (`template_creation.py`)

Demonstrates consensus template generation:
- Creating templates from site clusters
- Breaking templates into functional residue groups
- Analyzing template composition

**Usage:**
```bash
python examples/template_creation.py
```

**Output:**
- Template statistics
- Residue group analysis
- Functional residue identification

## Common Issues

### "Could not import Mcsa module"
- Make sure you're running examples from the CSA-3D root directory
- Check that all dependencies are installed: `pip install -r requirements.txt`

### "CSA-3D database not available"
- Set `CSA3D_DATA_DIR` environment variable
- Or download data from [M-CSA database](https://www.ebi.ac.uk/thornton-srv/m-csa/)

### Clustering or template errors
- Some features require additional dependencies (dendrogram-cut, scipy)
- Install dev dependencies: `pip install -r requirements-dev.txt`

## Extending the Examples

To create your own scripts:

1. Import the main Mcsa class:
   ```python
   from Mcsa import Mcsa
   db = Mcsa()
   ```

2. Build entries:
   ```python
   db.build(mcsa_id)
   entry = db.entries[mcsa_id]
   ```

3. Access sites and perform analysis:
   ```python
   for site in entry.pdbsites:
       print(site.id)  # Site identifier
       print(len(site))  # Number of residues
   ```

See the [API documentation](../docs/api_reference.md) for complete details.

## See Also

- [README.md](../README.md) - Installation and quick start
- [Architecture Guide](../docs/architecture.md) - System design
- [Contributing Guide](../CONTRIBUTING.md) - Development guidelines
