Building 3D Templates
=====================

Create consensus active site templates for matching.

Basic Template Creation
-----------------------

Build a template from conserved sites::

    template = entry.create_template(outdir="templates")
    
    # This writes a PDB file with representative coordinates of conserved residues.

The representative coordinates derive from the site with the highest
similarity to the average coordinates of the conserved residues.
The template nomenclature allows "fuzzy" matching of residues during template comparison.
This allows for some flexibility in matching residues that are not exactly identical but are chemically similar.
e.g. matching a Asp to a Glu due to their similar chemical properties.
You can view the exact definitions of fuzzy matching in the ``csa3d/residue_definitions.py`` module.

| The resulting template is suitable for 3D template matching with:

| - `Jess <https://github.com/iriziotis/jess>`_ (standalone tool)
| - `PyJess <https://github.com/althonos/pyjess>`_ (Python interface)

Custom Subsets
--------------

Build templates from specific sites::

    # Use only reference and homologous structures in a cluster
    template = entry.create_template(subset=[site.id for site in cluster_sites])
    
    # Use only a specific PDB structure
    template = entry.create_template(pdb_id="1ABC")

Template Parameters
-------------------

Customize template behavior::

    template = entry.create_template(
        ca=False,                      # False: all atoms, True: C-alpha only
        outdir="templates",            # Output directory
        outfile="template_1.pdb",      # Custom filename
        subset=None,                   # Site IDs to include (None = all)
        no_write=False,                # If True, don't write file
    )

Options:

- ``ca``: Use C-alpha atoms only (simpler, faster matching)
- ``outdir``: Where to write the template PDB file
- ``outfile``: Custom template filename (auto-named if None)
- ``subset``: Restrict to specific sites (list of site IDs)
- ``no_write``: Only compute template, don't write file

Template Breaking
------------------

Split large templates into smaller regions::

    template = entry.create_template()
    
    # Break into residue groups
    sub_templates = entry.break_template(
        template,
        method='clusters',
        n_residues=3,           # 3 residues per sub-template
        max_distance=6.0,       # Residues within 6 Å
    )
    
    # Write sub-templates
    for i, subtemp in enumerate(sub_templates):
        entry.write_template(subtemp, outdir="templates", 
                           outfile=f"sub_template_{i}.pdb")

Template Inspection
-------------------

Check template properties::

    print(f"Template ID: {template.id}")
    print(f"Number of residues: {len(template)}")
    print(f"Sequence: {template.sequence}")
    print(f"Resolution: {template.resolution}")
    
    for res in template.residues:
        print(f"{res.resname} {res.chain}{res.auth_resid}")

Writing Templates
-----------------

Write computed templates to disk::

    # Default location
    entry.write_template(template)
    
    # Custom location and name
    entry.write_template(
        template,
        outdir="/path/to/templates",
        outfile="my_template.pdb"
    )

Template Properties
-------------------

Templates inherit from PdbSite, so you can:

- Access residues: ``template.residues``
- Get distances: ``template.get_distances()``
- Check conservation: ``template.is_conserved``
- Examine ligands: ``template.ligands``

Integration with Jess
---------------------

Use CSA-3D templates with the Jess template-matching program::

    # CSA-3D: Build template
    template = entry.create_template(outdir="templates")
    
    # Jess/PyJess: Search PDB structures for similar arrangements (see relevant documentation)

See Also
--------

- :meth:`csa3d.Entry.create_template` for template creation API
- :meth:`csa3d.Entry.break_template` for splitting large templates
- :meth:`csa3d.Entry.write_template` for writing templates
- `Jess <https://github.com/iriziotis/jess>`_ for template matching
- `PyJess <https://github.com/althonos/pyjess>`_ for Python interface
