Comparing Catalytic Sites
==========================

Structural alignment and superposition of sites.

Basic Superposition
-------------------

Align two sites using the Superimposer algorithm::

    sites = list(entry.get_pdbsites(sane_only=True))
    site1 = sites[0]
    site2 = sites[1]
    
    # Superimpose site2 onto site1
    rot, tran, rms, rms_all = site1.fit(site2, weighted=True)
    
    print(f"Weighted RMSD: {rms:.2f} Ų")
    print(f"All atoms RMSD: {rms_all:.2f} Ų")

Parameters
----------

Tune the fitting algorithm::

    # Fast, robust fitting with outlier rejection
    rot, tran, rms, rms_all = site1.fit(
        site2,
        weighted=False,
        cycles=3,
        cutoff=3.0
    )
    
    # Precise fitting with Gaussian weighting
    rot, tran, rms, rms_all = site1.fit(
        site2,
        weighted=True,
        cycles=1,
        cutoff=999
    )
    
    # C-alpha only (for low-resolution or incomplete structures)
    rot, tran, rms, rms_all = site1.fit(site2, ca=True)

**Key parameters:**

- ``weighted``: Use Gaussian-weighted fitting (more precise) vs unweighted (faster)
- ``cycles``: Number of outlier rejection iterations
- ``cutoff``: Distance threshold for outlier rejection (Angstroms)
- ``ca``: Use only C-alpha atoms instead of all functional atoms
- ``transform``: If True, actually move site2's coordinates (default: False, compute only)

Return Values
-------------

The fit method returns::

    rot, tran, rms, rms_all = site1.fit(site2)
    
    # rot: 3×3 rotation matrix (NumPy array)
    # tran: 3-element translation vector (NumPy array)
    # rms: RMSD of atoms used in fitting (after outlier rejection)
    # rms_all: RMSD over all atoms

Apply the transformation manually::

    import numpy as np
    transformed = site2.coordinates @ rot + tran

Pairwise RMSD Matrix
--------------------

Compare all sites pairwise::

    sites = list(entry.get_pdbsites(sane_only=True))
    matrix = entry.rmsd_matrix(sites)
    
    # matrix is a pandas DataFrame
    print(matrix)

Per-Residue RMSD
----------------

See which residues differ most between sites::

    per_res_rms = site1.per_residue_rms(site2, rot=rot, tran=tran)
    
    for i, (res1, res2) in enumerate(zip(site1.residues, site2.residues)):
        if per_res_rms[i] > 2.0:  # Threshold
            print(f"Residue {res1.chain}{res1.auth_resid}: {per_res_rms[i]:.2f} Ų")

Geometric Comparison
--------------------

Compare local geometry without full alignment::

    # Inter-residue distances
    dist1 = site1.get_distances(kind='com')
    dist2 = site2.get_distances(kind='com')
    
    # Distance between specific residues
    d1 = site1.residues[0].get_distance(site1.residues[1])
    d2 = site2.residues[0].get_distance(site2.residues[1])

Functional Site Geometry
------------------------

Compare only catalytic atoms (excluding non-essential atoms)::

    func_site1 = site1.get_functional_site(ca=False)
    func_site2 = site2.get_functional_site(ca=False)
    
    rot, tran, rms, rms_all = func_site1.fit(func_site2, weighted=True)

Conservation Analysis
---------------------

Check if sites are structurally conserved::

    for site in entry.get_pdbsites():
        print(f"{site.id}:")
        print(f"  Conserved: {site.is_conserved}")
        print(f"  Conservative mutations: {site.is_conservative_mutation}")

Missing Atoms
-------------

Check for problems in structural quality::

    for site in entry.pdbsites:
        if site.has_missing_functional_atoms:
            print(f"{site.id} has missing functional atoms")
        if not site.is_sane:
            print(f"{site.id} failed sanity checks")

See Also
--------

- :meth:`csa3d.PdbSite.fit` for detailed fit() API
- :meth:`csa3d.Entry.rmsd_matrix` for pairwise comparisons
- :class:`csa3d.Superimposer` for underlying algorithm
- :doc:`clustering` for grouping similar sites
