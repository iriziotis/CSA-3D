# CSA-3D API Reference

Everything below is importable from the top-level package:

```python
import csa3d
from csa3d import Mcsa, Entry, PdbSite, PdbResidue, Het, UniSite, UniResidue
```

For a runnable tour of the main workflow, see
[`examples/basic_functionality.ipynb`](../examples/basic_functionality.ipynb).

---

## Data setup and configuration

### `csa3d.setup_data(data_dir=None, force=False) -> bool`

Download and prepare the reference data, then make `data_dir` the active data directory.

- `data_dir`: where to store the data. Defaults to `~/.csa3d/data`, which is found automatically in
  later sessions.
- `force`: re-download files that are already present. By default existing files are skipped, so
  re-running only fetches what is missing.

Returns `True` if every source was set up. Progress is reported through `logging`
(`logging.basicConfig(level=logging.INFO)` to see it). See the README for the list of sources.

### `csa3d.set_data_dir(data_dir: str) -> None`

Use an existing data directory for this session. Raises `ConfigError` if it does not exist.
Alternatively, set the environment variable `CSA3D_DATA_DIR`.

### `csa3d.get_data_dir() -> str`

The active data directory: the one set with `set_data_dir()` / `setup_data()`, else `$CSA3D_DATA_DIR`,
else `~/.csa3d/data`. Raises `ConfigError` if none exists.

### `csa3d.validate_config() -> bool`

Check that the core data files exist in the active data directory; raises `ConfigError` listing any
that are missing.

### `csa3d.datasets`

Lazily loaded lookup tables (SIFTS mappings, ion/cofactor sets, resolutions, ...), read on first access
from the active data directory and cached, e.g. `csa3d.datasets.PDB2EC[("1b73", "A")]`.

---

## `Mcsa`: the M-CSA database

```python
db = Mcsa()
```

Loads the M-CSA catalytic residue annotations and entry metadata from the data directory.

| Attribute | Description |
|---|---|
| `entries` | `{mcsa_id: Entry}` for every entry built so far |
| `catalytic_residue_info` | Raw M-CSA residue records, `{mcsa_id: [records]}`; its keys are all available entry IDs |
| `mcsa_entry_info` | Raw M-CSA entry metadata, `{mcsa_id: dict}` |

### `build(entry_ids=None, annotate=True, redundancy_cutoff=None, resolution_cutoff=None, verbose=False, no_sites=False) -> bool`

Build `Entry` objects with their PDB and UniProt catalytic sites. Assembly structures that are not yet in
the data directory are downloaded from PDBe.

- `entry_ids`: an M-CSA ID or list of IDs. `None` builds all entries (slow; downloads thousands of structures).
- `annotate`: find ligands near each site and read structure metadata (resolution, organism, ...).
- `redundancy_cutoff`: drop sites within this RMSD (Å) of another site in the same structure. With
  `None`, only identical sites are dropped.
- `resolution_cutoff`: skip non-reference structures whose resolution is worse than this (Å) or unknown.
- `verbose`: log each site as it is built.
- `no_sites`: only build residues, not sites.

Returns `False` if an ID is unknown or the entry has several reference structures (not supported).

### `update_pdbs(entry_ids=None) -> None`

Download the assembly structures for the given entries (default: all) ahead of time.

---

## `Entry`: one M-CSA entry

| Attribute | Description |
|---|---|
| `mcsa_id` | M-CSA identifier |
| `info` | M-CSA entry metadata: `enzyme_name`, `all_ecs`, `reaction`, `reference_uniprot_id`, ... |
| `pdbsites` | All `PdbSite`s (reference site and homologues) |
| `unisites` | All `UniSite`s |

| Method | Description |
|---|---|
| `get_pdbsites(pdb_id=None, sane_only=False)` | Iterate over PDB sites, optionally for one PDB ID and/or only sites without missing atoms or other problems |
| `get_pdbsite(id)`, `get_unisite(id)` | Look up a site by its `id` |
| `rmsd_matrix(sitelist, ca=False)` | Pairwise RMSD `pandas.DataFrame` (sites with missing functional atoms are skipped unless `ca=True`) |
| `clustering(matrix, height=None, get_linkage_matrix=False)` | Average-linkage hierarchical clustering; returns `{cluster_id: [site ids]}`. Cuts at 70% of the tree height by default |
| `clustering_bayesian(matrix, k_max=5, l=1.0, plot_outfile=None)` | Clustering with an automatic, Bayesian choice of the cut |
| `create_template(ca=False, outdir=None, outfile=None, subset=None, ...)` | Average the conserved sites into a consensus `PdbSite` and write it as a Jess template (unless `no_write=True`) |
| `write_template(template, ...)` | Write a template file for [Jess](https://github.com/iriziotis/jess) / [PyJess](https://github.com/althonos/pyjess) |
| `break_template(template, method='clusters', n_residues=3, max_distance=6.0)` | Split a large template into smaller residue groups |
| `update_mapped_sites()` | Link `UniSite`s and `PdbSite`s using SIFTS |

---

## `PdbSite`: a catalytic site in a PDB structure

| Attribute / property | Description |
|---|---|
| `id` | Unique site ID: PDB ID plus the chain and residue number of each catalytic residue |
| `pdb_id`, `uniprot_id`, `ec` | Identifiers (UniProt and EC via SIFTS) |
| `residues` | `PdbResidue`s, in the same order as the reference site (gaps where a residue has no equivalent) |
| `ligands` | `Het`s within 5 Å of the catalytic residues (when built with `annotate=True`) |
| `is_reference`, `is_conserved`, `is_conservative_mutation` | Relationship to the entry's reference site |
| `is_sane` | Whether the site passed sanity checks (see `Entry.get_pdbsites(sane_only=True)`) |
| `has_missing_functional_atoms` | Whether functional atoms are missing from the structure |
| `resolution`, `experimental_method`, `organism_name`, `title`, `enzyme`, `assembly_id` | Structure metadata (read from the mmCIF file) |
| `sequence`, `size` | One-letter pseudo-sequence and residue count |

### `fit(other, weighted=False, cycles=1, cutoff=999, transform=False, ca=False, ...) -> (rot, tran, rms, rms_all)`

Superimpose `other` (mobile) onto this site (fixed), finding the best atom correspondence within each
residue (accounting for symmetric atoms and conservative mutations). With `cycles > 1`, atoms further
than `cutoff` Å are excluded and the fit is repeated. Returns the rotation matrix, translation vector,
RMSD excluding outliers (`rms`) and RMSD over all atoms (`rms_all`). `transform=True` also moves
`other`'s coordinates.

### Other methods

| Method | Description |
|---|---|
| `per_residue_rms(other, rot=None, tran=None)` | RMSD of each residue after superposition |
| `get_functional_site(ca=False)` | Copy containing only the functional atoms of each residue |
| `write_pdb(outdir=None, outfile=None, write_hets=False, func_atoms_only=False)` | Write the site (optionally with ligands) as a PDB file |
| `get_residues()`, `get_gaps()` | Non-gap residues; gap residues |
| `get_distances(kind='com')` | Pairwise inter-residue distances |

---

## `PdbResidue`: a catalytic residue with coordinates

Key attributes: `resname`, `chain`, `resid`, `auth_resid`, `funclocs` (functional location, e.g.
side chain or main chain), `is_reference`, `reference_residue`, `structure` (the Biopython residue).
Properties: `is_gap`, `is_conserved`, `is_conservative_mutation`, `is_standard`, `has_main_chain_function`.

## `Het`: a ligand near a catalytic site

Key attributes: `resname`, `chain`, `resid`, `similarity` (fingerprint similarity to the enzyme's
reaction components, or `None`), `centrality`.
`type` classifies the component as artefact, polymer, cofactor, metal, ion or ligand; see also
`is_metal`, `is_ion`, `is_cofactor`, `is_artefact`, `is_polymer`.

---

## `UniSite` and `UniResidue`: sites at sequence level

`UniSite` mirrors `PdbSite` for UniProt sequences: `id` / `uniprot_id`, `residues`, `sequence`,
`is_reference`, `is_conserved`, `is_conservative_mutation`, `mapped_pdbsites` (after
`Entry.update_mapped_sites()`), `get_residues()`, `get_gaps()`.

`UniResidue` holds `resname`, `resid`, `funclocs` and `reference_residue`, with `is_gap`,
`is_conserved` and `is_conservative_mutation`.

---

## See also

- [README](../README.md): installation, data setup and usage
- [Architecture guide](./architecture.md)
- [Examples](../examples/README.md)
