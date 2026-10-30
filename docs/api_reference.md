# CSA-3D API Reference

Complete API documentation for CSA-3D classes and methods.

---

## Module: `csa3d`

Top-level package initialization and utilities.

### Function: `setup_data(data_dir: Optional[str] = None) -> bool`

**Description**: Automatically download and setup CSA-3D data files.

Downloads all required data from official sources and creates proper directory structure. This is the recommended way to initialize CSA-3D for first-time use.

**Parameters:**
- `data_dir` (str, optional): Directory for data files. Defaults to `~/.csa3d/data`

**Returns:** bool - True if setup successful, False if any downloads failed

**Examples:**

```python
from csa3d import setup_data

# Use default location
setup_data()

# Use custom location
setup_data("/custom/data/path")

# With error handling
if setup_data():
    print("✓ Data setup complete!")
else:
    print("⚠ Some data sources failed - check logs")
```

**Downloads:**
- SIFTS: PDB-UniProt mappings from EBI
- RHEA: Reaction component data
- PDBe: Ion, cofactor, and resolution data
- M-CSA: Catalytic site data

**Notes:**
- ~50 GB total download size
- One-time setup required
- Sets `CSA3D_DATA_DIR` environment variable automatically
- See `csa3d.data_setup` module for advanced usage

---

### Function: `set_data_dir(data_dir: str) -> None`

**Description**: Set the CSA-3D data directory manually.

Use this if you've already downloaded data or prefer manual setup.

**Parameters:**
- `data_dir` (str): Path to directory containing CSA-3D data files

**Returns:** None

**Raises:** ConfigError if directory does not exist

**Example:**

```python
from csa3d import set_data_dir

set_data_dir("/path/to/csa3d_data")
```

---

### Function: `get_data_dir() -> str`

**Description**: Get the current CSA-3D data directory.

**Returns:** str - Path to data directory

**Example:**

```python
from csa3d import get_data_dir

data_dir = get_data_dir()
print(f"Data directory: {data_dir}")
```

---

## Module: `Mcsa`

**Description**: Top-level database class for M-CSA entries.

Contains Entry objects and provides methods to load and build entries from M-CSA JSON data. Manages catalytic sites from PDB structures and UniProt sequences.

### Class: `Mcsa`

Main database class for accessing and building M-CSA entries.

**Key Methods:**

#### `__init__()`
Initialize the Mcsa database.

**Parameters:**
- None

**Returns:** None

---

#### `build(mcsa_id: int, annotate: bool = False, redundancy_cutoff: float = 0.95, resolution_cutoff: float = 2.5, verbose: bool = False) -> bool`

Build a specific M-CSA entry.

**Parameters:**
- `mcsa_id` (int): M-CSA entry identifier
- `annotate` (bool): Whether to include annotations
- `redundancy_cutoff` (float): Sequence redundancy threshold
- `resolution_cutoff` (float): PDB resolution threshold (Å)
- `verbose` (bool): Print progress information

**Returns:** bool - True if successfully built, False otherwise

---

#### `add(entry: Entry) -> bool`

Add an Entry to the database.

**Parameters:**
- `entry` (Entry): Entry object to add

**Returns:** bool - True if successfully added

---

#### `update_pdbs(pdbs: List[str]) -> None`

Update PDB structures in the database.

**Parameters:**
- `pdbs` (List[str]): List of PDB IDs to update

**Returns:** None

---

## Module: `Entry`

**Description**: M-CSA entry representation and analysis.

Provides the Entry class for manipulating enzyme catalytic sites from both PDB structures and UniProt sequences. Includes RMSD calculation, clustering, and template generation methods.

### Class: `Entry`

Represents an M-CSA entry with PDB and UniProt sites.

**Attributes:**
- `mcsa_id` (int): M-CSA identifier
- `e_number` (str): EC number
- `pdbsites` (List[PdbSite]): Catalytic sites from PDB structures
- `unisites` (List[UniSite]): Catalytic sites from UniProt sequences

**Key Methods:**

#### `__init__(mcsa_id: int, e_number: str) -> None`

Initialize an Entry.

**Parameters:**
- `mcsa_id` (int): M-CSA entry identifier
- `e_number` (str): EC enzyme number

**Returns:** None

---

#### `add(site: Union[PdbSite, UniSite]) -> None`

Add a catalytic site to the entry.

**Parameters:**
- `site` (Union[PdbSite, UniSite]): Site to add

**Returns:** None

---

#### `rmsd_matrix(sitelist: List[PdbSite]) -> np.ndarray`

Calculate pairwise RMSD matrix for sites.

**Parameters:**
- `sitelist` (List[PdbSite]): Sites to compare

**Returns:** np.ndarray - 2D RMSD matrix

---

#### `clustering(matrix: np.ndarray) -> Dict[int, List[int]]`

Perform hierarchical clustering on RMSD matrix.

**Parameters:**
- `matrix` (np.ndarray): RMSD matrix

**Returns:** Dict[int, List[int]] - Cluster ID mapped to site indices

---

#### `clustering_bayesian(matrix: np.ndarray, plot_outfile: Optional[str] = None) -> Dict[int, List[int]]`

Perform Bayesian hierarchical clustering.

**Parameters:**
- `matrix` (np.ndarray): RMSD matrix
- `plot_outfile` (Optional[str]): Path to save dendrogram

**Returns:** Dict[int, List[int]] - Cluster assignments

---

#### `create_template(ca: bool = True, outdir: Optional[str] = None, outfile: Optional[str] = None, subset: Optional[List[int]] = None, cluster_no: int = 0, atoms: str = "CA", no_write: bool = False, all_fuzzy: bool = False, no_alt: bool = False) -> PdbSite`

Create a consensus template from catalytic sites.

**Parameters:**
- `ca` (bool): Use CA atoms only
- `outdir` (Optional[str]): Output directory
- `outfile` (Optional[str]): Output filename
- `subset` (Optional[List[int]]): Site indices to use
- `cluster_no` (int): Cluster number
- `atoms` (str): Atom selection ("CA", "CB", etc.)
- `no_write` (bool): Don't write output file
- `all_fuzzy` (bool): Include fuzzy residues
- `no_alt` (bool): Skip alternate conformations

**Returns:** PdbSite - Template structure

---

#### `break_template(template: PdbSite, n_residues: int = 7, permutations: bool = False, max_distance: float = 8.0) -> Dict[int, Tuple[List[str], List[str]]]`

Divide template into functional residue groups.

**Parameters:**
- `template` (PdbSite): Template structure
- `n_residues` (int): Residues per group
- `permutations` (bool): Generate permutations
- `max_distance` (float): Distance threshold (Å)

**Returns:** Dict - Residue group mappings

---

## Module: `PdbSite`

**Description**: Protein Data Bank (PDB) structure site representation.

Provides the PdbSite class for managing catalytic sites from PDB structures, including residue extraction, RMSD calculation, and structure fitting.

### Class: `PdbSite`

Represents a catalytic site from a PDB structure.

**Attributes:**
- `id` (str): Site identifier
- `pdb_id` (str): PDB structure ID
- `residues` (List[PdbResidue]): Catalytic residues

**Key Methods:**

#### `__init__(pdb_id: str, site_id: str) -> None`

Initialize a PdbSite.

**Parameters:**
- `pdb_id` (str): PDB structure identifier
- `site_id` (str): Site identifier

**Returns:** None

---

#### `add(residue: PdbResidue) -> None`

Add a residue to the site.

**Parameters:**
- `residue` (PdbResidue): Residue to add

**Returns:** None

---

#### `fit(other: PdbSite, weighted: bool = False, ca: bool = True, scaling_factor: Optional[float] = None, transform: bool = False) -> Tuple[np.ndarray, np.ndarray, float, float]`

Fit this site to another structure.

**Parameters:**
- `other` (PdbSite): Target site
- `weighted` (bool): Use weighted fitting
- `ca` (bool): Use CA atoms only
- `scaling_factor` (Optional[float]): Scaling for coordinates
- `transform` (bool): Apply transformation

**Returns:** Tuple[rotation_matrix, translation_vector, rmsd, rmsd_all]

---

#### `write_pdb(outdir: str, outfile: Optional[str] = None, func_atoms_only: bool = False, write_hets: bool = True) -> None`

Write structure to PDB file.

**Parameters:**
- `outdir` (str): Output directory
- `outfile` (Optional[str]): Output filename
- `func_atoms_only` (bool): Write functional atoms only
- `write_hets` (bool): Include heteroatoms

**Returns:** None

---

## Module: `UniSite`

**Description**: UniProt catalytic site representation for M-CSA.

### Class: `UniSite`

Represents a catalytic site sequence from UniProt.

**Attributes:**
- `uniprot_id` (str): UniProt identifier
- `residues` (List[UniResidue]): Sequence residues

**Key Methods:**

#### `__init__(uniprot_id: str) -> None`

Initialize a UniSite.

**Parameters:**
- `uniprot_id` (str): UniProt identifier

**Returns:** None

---

#### `add(residue: UniResidue) -> None`

Add a residue to the site.

**Parameters:**
- `residue` (UniResidue): Residue to add

**Returns:** None

---

## Module: `PdbResidue`

**Description**: PDB catalytic residue representation for M-CSA.

### Class: `PdbResidue`

Represents a residue from a PDB structure.

**Attributes:**
- `residue_name` (str): Three-letter residue code
- `residue_number` (int): Position in structure
- `chain` (str): Chain identifier
- `atoms` (Dict): Atom coordinates

**Key Methods:**

#### `__init__(residue_name: str, residue_number: int, chain: str, ins_code: str = "") -> None`

Initialize a PdbResidue.

**Parameters:**
- `residue_name` (str): Three-letter code
- `residue_number` (int): Sequence number
- `chain` (str): Chain ID
- `ins_code` (str): Insertion code

**Returns:** None

---

## Module: `UniResidue`

**Description**: UniProt sequence residue representation for M-CSA catalytic sites.

### Class: `UniResidue`

Represents a residue from a UniProt sequence.

**Attributes:**
- `residue_name` (str): Single or three-letter code
- `position` (int): Position in sequence
- `is_conserved` (bool): Conservation status
- `function` (str): Functional role

**Key Methods:**

#### `__init__(residue_name: str, position: int) -> None`

Initialize a UniResidue.

**Parameters:**
- `residue_name` (str): Amino acid identifier
- `position` (int): Position in sequence

**Returns:** None

---

## Configuration Module

### Class: `DataPaths`

Manages data file paths for CSA-3D.

**Methods:**

#### `set_data_dir(data_dir: str) -> None`

Set the data directory explicitly.

**Parameters:**
- `data_dir` (str): Path to data directory

**Returns:** None

---

#### `get_data_dir() -> Path`

Get the configured data directory.

**Returns:** Path - Data directory path

---

## Environment Variables

### CSA3D_DATA_DIR

Path to directory containing M-CSA data files.

**Default locations:**
- `~/.csa3d/data`
- `/nfs/research/thornton/...` (EBI specific)
- User-specified via environment variable

**Example:**
```bash
export CSA3D_DATA_DIR=/path/to/data
```

---

## Common Data Types

### Dictionary Types

- **Clusters**: `Dict[int, List[int]]` - Maps cluster ID to site indices
- **Templates**: `Dict[int, Tuple[List[str], List[str]]]` - Residue groups
- **Mappings**: `Dict[str, str]` or similar - Data lookups

### Array Types

- **RMSD Matrix**: `np.ndarray` - 2D array of pairwise distances
- **Coordinates**: `np.ndarray` - 3D atom positions
- **Transformation**: `Tuple[np.ndarray, np.ndarray]` - Rotation and translation

---

## See Also

- [README.md](../README.md) - Installation and quick start
- [Architecture Guide](./architecture.md) - System design
- [Examples](../examples/README.md) - Usage examples
- [Contributing Guide](../CONTRIBUTING.md) - Development guidelines
