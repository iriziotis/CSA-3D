Data Setup and Configuration
==============================

Managing reference data and custom data directories.

Initial Setup
-------------

CSA-3D requires reference data (~400 MB) covering:

- M-CSA catalytic residue annotations
- PDB structure metadata and assemblies
- SIFTS UniProt-PDB mappings
- RHEA and KEGG reaction component data
- Ion and cofactor databases

Download everything with::

    import csa3d
    csa3d.setup_data()

This creates ``~/.csa3d/data`` (or ``$CSA3D_DATA_DIR`` if set).

Custom Data Directory
---------------------

Store data in a custom location::

    csa3d.setup_data("/path/to/data")

Or set the location for this session::

    csa3d.set_data_dir("/path/to/data")

And for all future sessions, set the environment variable::

    export CSA3D_DATA_DIR=/path/to/data

Checking Your Setup
-------------------

Verify data is complete::

    csa3d.validate_config()

This checks that all required data files exist and raises ``ConfigError``
if anything is missing.

Structure Downloads
-------------------

Structure files (mmCIF) are downloaded automatically when needed during
``Mcsa.build()``. Pre-download to work offline::

    db = Mcsa()
    db.update_pdbs(entry_ids=[1, 2, 3])

Or download all PDB structures for all M-CSA entries (slow)::

    db.update_pdbs()

Bundled Data
------------

Some reference data ships with CSA-3D and is copied to the data directory:

- ``mcsa/mcsa_cofactors.csv`` — M-CSA cofactor annotations
- ``pdbe/pdb_ions.csv`` — Common PDB ions
- ``pdbe/pdb_cofactors.csv`` — Common PDB cofactors
- ``pdbe/crystallization_hets.csv`` — Crystallographic artifacts
- ``pdbe/metals.csv`` — Metal and metal-containing compounds

These snapshots are used if remote sources are unavailable.

Updating Data
-------------

Files are downloaded on-demand and cached. To refresh everything::

    csa3d.setup_data(force=True)

To only update missing files::

    csa3d.setup_data()

Data Directory Contents
-----------------------

After setup, your data directory contains::

    data/
    ├── mcsa/                          # M-CSA database
    │   ├── catalytic_residues_references.json
    │   ├── catalytic_residues_homologues.json
    │   ├── entry_info/
    │   └── mcsa_cofactors.csv
    ├── sifts/                         # UniProt-PDB mappings
    │   ├── pdb_chain_taxonomy.csv
    │   └── *.csv (one per PDB chain)
    ├── pdbe/                          # PDB metadata
    │   ├── resolu.idx                 # Structure resolutions
    │   ├── pdb_ions.csv
    │   ├── pdb_cofactors.csv
    │   ├── crystallization_hets.csv
    │   ├── metals.csv
    │   └── assembly_cif/cif/          # Structures (downloaded on demand)
    ├── rhea/                          # RHEA reaction data
    │   └── *.tsv
    └── kegg/                          # KEGG reaction data
        └── *.csv

Environment Variables
---------------------

- ``CSA3D_DATA_DIR`` — Path to data directory (auto-discovered in ``~/.csa3d/data``)

See Also
--------

- :func:`csa3d.setup_data`
- :func:`csa3d.set_data_dir`
- :func:`csa3d.get_data_dir`
- :func:`csa3d.validate_config`
