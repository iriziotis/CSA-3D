Quick Start
===========

Get up and running with CSA-3D in 5 minutes.

1. Installation
---------------

::

    pip install csa3d

2. Download Reference Data
---------------------------

One-off setup (~400 MB download)::

    import csa3d
    csa3d.setup_data()

This downloads the M-CSA database, PDB structures, SIFTS mappings, and other
reference data to ``~/.csa3d/data``. In later sessions, CSA-3D automatically
finds data in this directory.

To use a custom location::

    csa3d.setup_data("/path/to/data")

3. Load the Database
--------------------

::

    from csa3d import Mcsa
    
    db = Mcsa()
    db.build(entry_ids=1)
    entry = db.entries[1]

This loads M-CSA entry 1 (Glutamate Racemases). The ``build()`` method constructs
``Entry`` objects with PDB structures and UniProt sequences automatically
downloaded on demand.

4. Explore Catalytic Sites
---------------------------

::

    # Get entry metadata
    print(entry.info["enzyme_name"])      # e.g., "Glutamate Racemase"
    print(entry.info["all_ecs"])          # e.g., ["5.1.1.3"]
    
    # Iterate over PDB sites
    for site in entry.get_pdbsites(sane_only=True):
        print(f"{site.id}: {site.pdb_id} resolution={site.resolution}")
        for res in site.residues:
            print(f"  {res.resname} {res.chain}{res.auth_resid}")

5. Compare Sites
----------------

Superimpose two catalytic sites::

    sites = list(entry.get_pdbsites(sane_only=True))
    rot, tran, rms, rms_all = sites[0].fit(sites[1], weighted=True)
    print(f"RMSD: {rms:.2f} Å (all atoms: {rms_all:.2f} Å)")

The ``fit()`` method uses the paper's Superimposer algorithm to find optimal
structural alignment.

6. Cluster by Geometry
----------------------------

Group sites by structural similarity::

    matrix = entry.rmsd_matrix(sites)
    clusters = entry.clustering(matrix)
    
    for cluster_id, site_ids in clusters.items():
        print(f"Cluster {cluster_id}: {len(site_ids)} sites")

7. Build a 3D Template
-----------------------

Create a consensus template for the catalytic site::

    template = entry.create_template(outdir="templates")

This writes a PDB file suitable for template-matching with
`Jess <https://github.com/iriziotis/jess>`_ or
`PyJess <https://github.com/althonos/pyjess>`_.

Next Steps
----------

- See :doc:`usage/data_setup` for data management options
- See :doc:`usage/compare_sites` for detailed comparison workflows
- See :doc:`api/core` for complete API documentation
- Check out the `examples <https://github.com/iriziotis/CSA-3D/tree/master/examples>`_ directory

Troubleshooting
---------------

**ImportError: No module named csa3d**
    Make sure CSA-3D is installed: ``pip install csa3d``

**ConfigError: Data directory not found**
    Run ``csa3d.setup_data()`` first, or set ``CSA3D_DATA_DIR`` environment
    variable.

**Cannot download PDB structure**
    Check your internet connection. The first build of an entry with many
    structures may take a while. Use ``db.update_pdbs()`` to pre-download.

For more help, see `GitHub Issues <https://github.com/iriziotis/CSA-3D/issues>`_.
