# CSA-3D: 3D Analysis of Enzyme Catalytic Sites

CSA-3D extracts, compares and manipulates enzyme catalytic sites in three dimensions,
using the catalytic residue annotations of the
[Mechanism and Catalytic Site Atlas (M-CSA)](https://www.ebi.ac.uk/thornton-srv/m-csa/).

With CSA-3D you can:

- **Extract** the catalytic site of every PDB structure (and UniProt sequence) homologous to an M-CSA entry
- **Compare** catalytic sites by structural superposition and RMSD
- **Cluster** sites by geometry to explore conformational and functional variation
- **Build 3D templates** that summarise an active site, for use with the template-matching software
  [Jess](https://github.com/iriziotis/jess) / [PyJess](https://github.com/althonos/pyjess)
- **Annotate** ligands, cofactors and metal ions found in catalytic sites

See the [CSA-3D paper](https://doi.org/10.1016/j.jmb.2022.167517) for the method and applications.

## Installation

Requires Python 3.8 or later.

```bash
pip install csa3d
```

or, from source:

```bash
git clone https://github.com/iriziotis/CSA-3D.git
cd CSA-3D
pip install -e .
```

## Quick start

```python
import csa3d

# 1. One-off: download the reference data (~400 MB) to ~/.csa3d/data
csa3d.setup_data()

# 2. Load M-CSA and build an entry
db = csa3d.Mcsa()
db.build(entry_ids=1)
entry = db.entries[1]
print(entry.info["enzyme_name"], entry.info["all_ecs"])

# 3. Explore its catalytic sites
for site in entry.get_pdbsites(sane_only=True):
    print(site.id, site.pdb_id, site.resolution, [res.resname for res in site.residues])
```

In later sessions you can skip step 1: `csa3d.Mcsa()` finds data in `~/.csa3d/data` automatically.
For a complete, runnable walkthrough see [`examples/basic_functionality.ipynb`](examples/basic_functionality.ipynb).

## Data

### Setting up

`csa3d.setup_data()` downloads and prepares everything CSA-3D needs, then makes that directory the
active one. It is safe to re-run: files already present are skipped, so a re-run only retries what
failed. Pass `force=True` to refresh everything.

```python
csa3d.setup_data()                    # default location: ~/.csa3d/data
csa3d.setup_data("/data/csa3d")       # or any location you choose
```

To see progress, enable logging first: `import logging; logging.basicConfig(level=logging.INFO)`.
From the command line: `python -m csa3d.data_setup /data/csa3d`.

### Using a custom location

The default location is found automatically. For any other location, tell CSA-3D where the data lives
at the start of each session, in one of two ways:

```python
csa3d.set_data_dir("/data/csa3d")
```

```bash
export CSA3D_DATA_DIR=/data/csa3d
```

### Protein structures

Biological assembly structures (mmCIF) are **downloaded from PDBe on demand** when `Mcsa.build()`
needs them, and cached in `<data dir>/pdbe/assembly_cif/cif/`. The first build of an entry with many
homologous structures therefore takes longer. To fetch structures up front (e.g. before working
offline), use `db.update_pdbs(entry_ids=[1, 2, 3])`, or `db.update_pdbs()` for all of M-CSA.

### What gets downloaded

| File (under the data directory) | Source |
|---|---|
| `mcsa/catalytic_residues_homologues.json` | M-CSA API: catalytic residues in homologous structures (~325 MB) |
| `mcsa/catalytic_residues_references.json` | M-CSA API: catalytic residues of reference structures |
| `mcsa/entry_info/mcsa_entry_info.*.json` | M-CSA API: entry metadata (enzyme name, EC, reaction) |
| `sifts/*.csv` | [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/): PDB–UniProt–EC mappings |
| `pdbe/resolu.idx` | wwPDB: resolution of every PDB entry |
| `rhea/rhea_reaction_components.csv` | Compiled from the latest [RHEA](https://www.rhea-db.org/) release |
| `kegg/kegg_reaction_components.csv` | Compiled from the [KEGG](https://www.kegg.jp/) REST API |
| `pdbe/assembly_cif/cif/*.cif` | PDBe, fetched on demand by `Mcsa.build()` |

Some datasets are no longer published, so snapshots **ship with the package** and are copied into the
data directory: `mcsa/mcsa_cofactors.csv`, `pdbe/pdb_ions.csv`, `pdbe/pdb_cofactors.csv`,
`pdbe/crystallization_hets.csv` and `pdbe/metals.csv`. RHEA and KEGG also fall back to bundled
snapshots if their servers cannot be reached.

Ligand similarity scoring (`Het.similarity`) additionally uses small-molecule structure files in
`pdbe/het_sdf/`, `chebi/sdf/` and `kegg/sdf/`, which `setup_data()` does not download. Without them,
`similarity` is `None`; everything else works.

## Usage

The examples below continue from the quick start (`entry = db.entries[1]`).

### Inspect a catalytic site

Each entry has one **reference site** (the structure annotated by M-CSA curators) and many homologous
sites mapped onto it.

```python
reference = next(site for site in entry.get_pdbsites() if site.is_reference)
print(reference.pdb_id, reference.uniprot_id, reference.ec, reference.resolution)

for res in reference.residues:
    print(res.resname, res.chain, res.auth_resid, res.funclocs)

for het in reference.ligands:          # ligands within 5 Å of the catalytic residues
    print(het.resname, het.type)
```

Use `entry.get_pdbsites(sane_only=True)` to skip sites with missing atoms or other problems.

### Compare two sites

```python
sites = list(entry.get_pdbsites(sane_only=True))
rot, tran, rms, rms_all = sites[0].fit(sites[1], weighted=True)
print(f"RMSD {rms:.2f} Å (all atoms: {rms_all:.2f} Å)")
```

`fit()` superimposes the second (mobile) site onto the first (fixed) one. `rms` excludes atoms rejected
as outliers during iterative fitting (`cycles`, `cutoff`); `rms_all` includes every atom.

### Cluster sites by geometry

```python
matrix = entry.rmsd_matrix(sites)      # pandas DataFrame of pairwise RMSDs
clusters = entry.clustering(matrix)    # {cluster_id: [site ids]}
for cluster_id, site_ids in clusters.items():
    print(cluster_id, len(site_ids))
```

`clustering()` uses average-linkage hierarchical clustering, cut at 70% of the tree height unless you
pass `height=`.

### Build a 3D template

```python
template = entry.create_template(outdir="templates")
```

This averages the conserved sites into a consensus geometry and writes a template file to
`templates/`. Templates are input for [Jess](https://github.com/iriziotis/jess), a 3D template-matching
program that searches protein structures for arrangements of atoms resembling the template, or for
its Python interface [PyJess](https://github.com/althonos/pyjess). Use `ca=True` for C-alpha templates
and `subset=[...]` to build from selected sites only (e.g. one cluster).

### Export structures

```python
reference.write_pdb(outdir="sites", write_hets=True)   # catalytic residues (and ligands) as PDB
```

### Sequence-level sites

```python
for unisite in entry.unisites:
    print(unisite.id, unisite.sequence, unisite.is_conserved)
```

### Build several entries

```python
db.build(entry_ids=[1, 2, 3], resolution_cutoff=2.5)
for entry in db:
    print(entry.mcsa_id, len(entry.pdbsites))
```

`db.build()` with no arguments builds all of M-CSA (~1000 entries); this downloads many thousands of
structures and takes a long time.

## Main classes

| Class | Represents | Key members |
|---|---|---|
| `Mcsa` | The M-CSA database | `build()`, `entries`, `update_pdbs()` |
| `Entry` | One M-CSA entry | `info`, `pdbsites`, `unisites`, `get_pdbsites()`, `rmsd_matrix()`, `clustering()`, `create_template()` |
| `PdbSite` | A catalytic site in a PDB structure | `residues`, `ligands`, `fit()`, `write_pdb()`, `is_reference`, `is_conserved` |
| `PdbResidue` | A catalytic residue with coordinates | `resname`, `chain`, `auth_resid`, `funclocs`, `structure` |
| `Het` | A ligand, cofactor or ion near a site | `resname`, `type`, `similarity` |
| `UniSite`, `UniResidue` | The same, at sequence level | `sequence`, `is_conserved` |

All are importable directly from `csa3d`, e.g. `from csa3d import Mcsa, PdbSite`.
See [docs/api_reference.md](docs/api_reference.md) for the full API.

## Testing

```bash
pip install -e ".[dev]"
pytest tests/
```

Running [`examples/basic_functionality.ipynb`](examples/basic_functionality.ipynb) end to end checks the
main workflow against real data.

## Citation

If you use CSA-3D in your research, please cite:

```bibtex
@article{RIZIOTIS2022167517,
    title = {Conformational Variation in Enzyme Catalysis: A Structural Study on Catalytic Residues},
    journal = {Journal of Molecular Biology},
    volume = {434},
    number = {7},
    pages = {167517},
    year = {2022},
    issn = {0022-2836},
    doi = {https://doi.org/10.1016/j.jmb.2022.167517},
    url = {https://www.sciencedirect.com/science/article/pii/S0022283622000912},
    author = {Ioannis G. Riziotis and António J.M. Ribeiro and Neera Borkakoti and Janet M. Thornton}
}
```

## Contributing

Contributions are welcome. Please open an issue or pull request on
[GitHub](https://github.com/iriziotis/CSA-3D); see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

Apache License 2.0; see [LICENSE](LICENSE).

## Contact

- **Author**: Ioannis Riziotis
- **Email**: riziotis@ebi.ac.uk
- **Bug reports**: [GitHub Issues](https://github.com/iriziotis/CSA-3D/issues)

CSA-3D builds on [Biopython](https://biopython.org/), [NumPy](https://numpy.org/),
[SciPy](https://scipy.org/), [pandas](https://pandas.pydata.org/), [scikit-learn](https://scikit-learn.org/),
[RDKit](https://www.rdkit.org/) and [rmsd](https://github.com/charnley/rmsd).
See [CHANGELOG.md](CHANGELOG.md) for version history.
