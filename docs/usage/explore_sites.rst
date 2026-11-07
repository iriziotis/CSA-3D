Exploring Catalytic Sites
===========================

Inspect and analyze individual catalytic sites.

Basic Site Properties
---------------------

After building an entry, access its sites::

    entry = db.entries[1]
    
    # Iterate over all sites
    for site in entry.pdbsites:
        print(f"Site {site.id}")
        print(f"  PDB ID: {site.pdb_id}")
        print(f"  Resolution: {site.resolution}")
        print(f"  Size: {len(site)} residues")

Filtering Sites
---------------

Get only high-quality structures::

    # Skip structures with missing atoms
    for site in entry.get_pdbsites(sane_only=True):
        print(site.id)
    
    # Get a specific PDB structure
    for site in entry.get_pdbsites(pdb_id="1E3E"):
        print(site)
    
    # Combine filters
    ref_site = next(s for s in entry.get_pdbsites(sane_only=True) if s.is_reference)

Residue Inspection
------------------

Examine catalytic residues::

    site = next(s for s in entry.pdbsites if s.is_reference)
    
    for res in site.residues:
        print(f"{res.resname} {res.chain}{res.auth_resid}")
        print(f"  Functional role: {res.funclocs}")
        print(f"  Is conserved: {res.is_conserved}")
        print(f"  Is reference: {res.is_reference}")
    
    # Get residue by ID
    res = site[("1E3E", "ALA", "A", 41, "L", 1)]

Ligands and Cofactors
---------------------

Find molecules near the active site::

    for het in site.ligands:
        print(f"{het.resname} (chain {het.chain})")
        print(f"  Type: {het.type}")
        print(f"  Similarity: {het.similarity}")
        print(f"  Centrality: {het.centrality}")
    
    # Filter by type
    metals = [h for h in site.ligands if h.is_metal]
    cofactors = [h for h in site.ligands if h.is_cofactor]
    artifacts = [h for h in site.ligands if h.is_artefact]

Sequence Level (UniProt)
------------------------

Examine sequence-level information::

    for unisite in entry.unisites:
        print(f"UniProt {unisite.uniprot_id}")
        print(f"  Sequence: {unisite.sequence}")
        print(f"  Is conserved: {unisite.is_conserved}")
        
        for res in unisite.residues:
            if not res.is_gap:
                print(f"    {res.resname} at position {res.resid}")

Structural Information
----------------------

Access structure-specific metadata::

    # Resolution and method
    print(site.resolution)
    print(site.experimental_method)
    print(site.organism_name)
    print(site.title)  # PDB entry title
    
    # Biological assembly
    print(site.assembly_id)
    
    # Associated enzyme information
    print(site.enzyme)  # If available
    print(site.ec)      # EC number (via SIFTS)

Writing Structures
------------------

Export sites as PDB files::

    # Write catalytic residues
    site.write_pdb(outdir="sites")
    
    # Include ligands
    site.write_pdb(outdir="sites", write_hets=True)
    
    # Only functional atoms (CA, CB, etc.)
    site.write_pdb(outdir="sites", func_atoms_only=True)

Distances and Geometry
----------------------

Calculate inter-residue distances::

    distances = site.get_distances(kind='com')  # Center of mass
    distances = site.get_distances(kind='cog')  # Center of geometry
    distances = site.get_distances(kind='ca')   # C-alpha
    
    # Distance between two residues
    res1 = site.residues[0]
    res2 = site.residues[1]
    dist = res1.get_distance(res2, kind='com')

See Also
--------

- :doc:`compare_sites` for structural alignment
- :class:`csa3d.PdbSite` for complete API
- :class:`csa3d.UniSite` for sequence-level sites
