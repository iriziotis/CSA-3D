CSA-3D Documentation
====================

**CSA-3D** extracts, compares and manipulates enzyme catalytic sites in three dimensions,
using the catalytic residue annotations of the
`Mechanism and Catalytic Site Atlas (M-CSA) <https://www.ebi.ac.uk/thornton-srv/m-csa/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   usage/data_setup
   usage/explore_sites
   usage/compare_sites
   usage/clustering
   usage/templates

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/core
   api/structures
   api/analysis
   api/utilities

.. toctree::
   :maxdepth: 2
   :caption: Reference

   architecture

Quick Start
-----------

Install CSA-3D::

    pip install csa3d

Download and prepare the reference data (one-off)::

    import csa3d
    csa3d.setup_data()

Load M-CSA and explore catalytic sites::

    db = csa3d.Mcsa()
    db.build(entry_ids=1)
    entry = db.entries[1]
    
    for site in entry.get_pdbsites(sane_only=True):
        print(site.id, site.pdb_id, site.resolution)

For a complete walkthrough, see :doc:`quickstart`.

Key Classes
-----------

- **Mcsa**: The M-CSA database with methods to build and query entries
- **Entry**: A single M-CSA catalytic mechanism entry
- **PdbSite**: A catalytic site from a PDB structure
- **UniSite**: A catalytic site from a UniProt sequence
- **Superimposer**: Structural superposition with weighted fitting

See the :doc:`API Reference <api/core>` for complete documentation.

Resources
---------

- `GitHub Repository <https://github.com/iriziotis/CSA-3D>`_
- `M-CSA Database <https://www.ebi.ac.uk/thornton-srv/m-csa/>`_
- `CSA-3D Paper <https://doi.org/10.1016/j.jmb.2022.167517>`_

Citation
--------

If you use CSA-3D in your research, please cite:

    Riziotis, I. G., Ribeiro, A. J. M., Borkakoti, N., & Thornton, J. M. (2022).
    Conformational Variation in Enzyme Catalysis: A Structural Study on Catalytic Residues.
    *Journal of Molecular Biology*, 434(7), 167517.
    https://doi.org/10.1016/j.jmb.2022.167517

License
-------

Apache License 2.0. See LICENSE for details.
