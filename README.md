# CSA-3D: 3D Analysis of Enzyme Catalytic Sites

A Python package for extracting, analyzing, and manipulating enzyme catalytic sites in three dimensions from the **Mechanism and Catalytic Site Atlas (M-CSA)** hosted at [EMBL-EBI](https://www.ebi.ac.uk/thornton-srv/m-csa/).

## What is CSA-3D?

CSA-3D enables computational biochemists and structural biologists to:

- **Extract** catalytic site geometries from PDB structures and UniProt sequences
- **Compare** active sites across homologous proteins using RMSD clustering
- **Create** structure templates for enzyme active sites with geometric constraints
- **Analyze** ligand proximity and chemical similarity in active sites
- **Visualize** and export catalytic site information in multiple formats

This package was developed to support research on enzyme mechanisms and is particularly useful for:
- Comparative structural analysis of enzyme families
- Active site template generation for molecular docking
- Catalytic site conservation studies
- Structure-guided annotation of enzyme function

**See also**: [CSA-3D paper](https://www.sciencedirect.com/science/article/pii/S0022283622000912)

## Installation

### Requirements

- Python 3.8 or higher
- ~50 GB of disk space for complete M-CSA dataset (can work with subsets)

### Installation from PyPI

```bash
pip install csa3d
```

### Installation from Source

```bash
git clone https://github.com/your-username/csa3d.git
cd csa3d
pip install -e ".[dev]"
```

### Development Installation

For development with testing and documentation tools:

```bash
pip install -e ".[dev,docs]"
```

## Quick Start

### 1. Automatic Data Setup (Recommended)

Download and setup M-CSA data automatically:

```python
from csa3d import setup_data

# Downloads ~50 GB of data to ~/.csa3d/data
# Set up specific location with:
# setup_data("/path/to/csa3d_data")
success = setup_data()

if success:
    print("✓ Data setup complete!")
```

Or via command line:
```bash
python -m csa3d.data_setup /path/to/csa3d_data
```

### 2. Manual Data Configuration

If you prefer to download data manually:

**Option A: Environment variable**
```bash
export CSA3D_DATA_DIR="/path/to/csa3d_data"
```

**Option B: Set in code**
```python
from csa3d import set_data_dir
set_data_dir("/path/to/csa3d_data")
```

### 3. Load and Inspect Active Sites

```python
from csa3d import Mcsa

# Initialize the database
db = Mcsa()

# Build a specific entry (e.g., M-CSA ID 1)
db.build(entry_ids=1, annotate=True)

# Iterate over entries
for entry in db:
    print(f"Entry {entry.mcsa_id}: {len(entry.pdbsites)} PDB sites")
    
    # Get information
    print(f"  EC Number: {entry.info['ec']}")
    print(f"  Enzyme: {entry.info['enzyme_name']}")
    
    # Inspect active sites
    for site in entry.get_pdbsites(sane_only=True):
        print(f"    PDB: {site.pdb_id}, {len(site.residues)} residues")
        print(f"    Resolution: {site.resolution} Å")
```

### 4. Clustering and Template Creation

```python
# Build multiple entries
db.build(entry_ids=[1, 2, 3], annotate=True)
entry = db.entries[1]

# Calculate RMSD between all active sites
sites = list(entry.get_pdbsites(sane_only=True))
rmsd_matrix = entry.rmsd_matrix(sites)

# Cluster sites by similarity
clusters = entry.clustering(rmsd_matrix, height=1.5)

print(f"Found {len(clusters)} clusters")
for cluster_id, site_ids in clusters.items():
    print(f"  Cluster {cluster_id}: {len(site_ids)} sites")

# Create a template from conserved sites
template = entry.create_template(outdir="./templates")
print(f"Template created: {template.id}")
```

## Data Requirements

CSA-3D requires data from several external sources:

### Core M-CSA Data
- `mcsa/catalytic_residues_homologues.json` - Homologous catalytic residues
- `mcsa/catalytic_residues_references.json` - Reference catalytic residues
- `mcsa/entry_info/mcsa_entry_info.*.json` - Entry metadata

### PDB Data
- `pdbe/assembly_cif/cif/` - PDB assembly structures (mmCIF format)
- `pdbe/pdb_ions.csv` - Ion classification
- `pdbe/pdb_cofactors.csv` - Cofactor information
- `pdbe/crystallization_hets.csv` - Crystallization artifacts

### Sequence/Mapping Data (SIFTS)
- `sifts/uniprot_pdb.csv` - UniProt to PDB mapping
- `sifts/pdb_chain_uniprot.csv` - PDB to UniProt mapping
- `sifts/pdb_chain_enzyme.csv` - PDB to EC number mapping

### Reaction Data
- `kegg/kegg_reaction_components.csv` - KEGG enzyme reactions
- `rhea/rhea_reaction_components.csv` - RHEA reactions
- `pdbe/resolu.idx` - PDB resolution index

### Molecule Data
- `kegg/sdf/` - KEGG small molecules
- `chebi/sdf/` - ChEBI compounds
- `pdbe/het_sdf/` - PDB heteroatom structures

To download all data:
```python
from csa3d.config import setup_data
setup_data("/path/to/csa3d_data")
```

Alternatively, download manually from:
- M-CSA: https://www.ebi.ac.uk/thornton-srv/m-csa/download/
- PDBe: https://www.ebi.ac.uk/pdbe/download/
- SIFTS: https://www.ebi.ac.uk/pdbe/docs/sifts/
- KEGG: https://www.kegg.jp/
- RHEA: https://www.rhea-db.org/

## Usage Examples

### Comparing Two Active Sites

```python
from csa3d import Mcsa

db = Mcsa()
db.build(entry_ids=1)
entry = db.entries[1]

sites = list(entry.get_pdbsites(sane_only=True))
if len(sites) >= 2:
    site1, site2 = sites[0], sites[1]
    
    # Structural alignment
    rotation, translation, rmsd_ca, rmsd_all = site1.fit(site2, weighted=True)
    print(f"RMSD (backbone): {rmsd_ca:.2f} Å")
    print(f"RMSD (all atoms): {rmsd_all:.2f} Å")
```

### Writing Templates for Docking

```python
# Create template with specific parameters
template = entry.create_template(
    ca=False,                    # Use all heavy atoms
    outdir="./templates",
    outfile="myenzyme.template.pdb",
    all_fuzzy=False,            # Use conservative fuzzy matching
    no_alt=False                # Include alternative residue options
)

# Template can be used with docking software like TESS/Jess
```

### Batch Processing

```python
from csa3d import Mcsa

db = Mcsa()

# Build all entries (takes time)
entry_ids = list(db.catalytic_residue_info.keys())
db.build(entry_ids=entry_ids[:10], verbose=True)  # First 10

# Analyze
for entry in db:
    n_sites = len(entry.pdbsites)
    n_clusters = len(entry.clustering(entry.rmsd_matrix(list(entry.get_pdbsites()))))
    print(f"{entry.mcsa_id}: {n_sites} sites, {n_clusters} clusters")
```

## Key Classes

### Mcsa
Main database class. Load and manage M-CSA entries.
```python
db = Mcsa()
db.build(entry_ids=1)  # Build specific entry
entry = db.entries[1]  # Access Entry object
```

### Entry
Represents one M-CSA entry containing PDB and UniProt active sites.
```python
entry.pdbsites           # List of PdbSite objects
entry.unisites           # List of UniSite objects
entry.create_template()  # Generate template
entry.clustering()       # RMSD-based clustering
```

### PdbSite
A catalytic site from a PDB structure.
```python
site.pdb_id              # PDB identifier
site.residues            # Catalytic residues
site.ligands             # Proximal ligands
site.fit(other_site)     # Structural alignment
```

### UniSite
A catalytic site from UniProt sequence.
```python
site.uniprot_id          # UniProt accession
site.residues            # Catalytic residues
```

## Configuration

### Environment Variables

- `CSA3D_DATA_DIR` - Path to M-CSA data directory (default: auto-detect)

### Configuration File

Create `~/.csa3d/config.yml`:
```yaml
data_dir: /path/to/csa3d_data
```

### In Code

```python
from csa3d import config

config.set_data_dir("/path/to/data")
config.validate_config()  # Check if valid
print(config.get_data_dir())
```

## Documentation

- [Architecture Guide](docs/architecture.md) - System design and data flow
- [API Reference](https://csa3d.readthedocs.io) - Full API documentation
- [Examples](examples/) - Jupyter notebooks and scripts

## Testing

Run tests with pytest:
```bash
pytest tests/
```

Run with coverage:
```bash
pytest --cov=csa3d tests/
```

## Performance Considerations

- **First run** is slow (~minutes) as PDB structures are parsed and cached
- **RMSD clustering** on large site sets can be memory-intensive
- **Template creation** involves fitting multiple structures (depends on cluster size)

Tips:
- Use `sane_only=True` to filter low-quality sites
- Process entries in batches
- Use `resolution_cutoff` to filter low-resolution structures
- Cache parsed entries with pickle for reuse

## Citation

If you use CSA-3D in your research, please cite:

```bibtex
@article{Riziotis2022,
    title={CSA-3D: {A} framework to manipulate enzyme catalytic sites in three dimensions},
    author={Riziotis, Ioannis and others},
    journal={Journal of Molecular Biology},
    volume={434},
    number={5},
    pages={167488},
    year={2022},
    issn={0022-2836},
    doi={10.1016/j.jmb.2022.167488},
    url={https://www.sciencedirect.com/science/article/pii/S0022283622000912}
}
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make changes and add tests
4. Run quality checks: `make lint`
5. Submit a pull request

## Issues & Support

- **Bug reports**: [GitHub Issues](https://github.com/your-username/csa3d/issues)
- **Questions**: Create a discussion or post an issue
- **Feature requests**: Open an issue with `[FEATURE]` prefix

## License

Apache License 2.0 - see LICENSE file for details

## Contact

- **Author**: Ioannis Riziotis
- **Email**: riziotis@ebi.ac.uk
- **Affiliation**: EMBL-EBI, Wellcome Trust Genome Campus, Cambridge, UK

## Acknowledgments

CSA-3D is built on top of several excellent packages:
- [Biopython](https://biopython.org/) - Biological computation
- [NumPy](https://numpy.org/) - Numerical computing
- [SciPy](https://scipy.org/) - Scientific computing
- [scikit-learn](https://scikit-learn.org/) - Machine learning
- [RDKit](https://www.rdkit.org/) - Cheminformatics

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for version history and breaking changes.

---

**Questions?** See [docs/faq.md](docs/faq.md) or open an issue.
