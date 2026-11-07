from pymol import cmd, stored
from pymol.constants_palette import palette_dict
import numpy as np

@cmd.extend
def color_catres(selection, palette=None):
    """ For each object in selection, assigns a distinct color to each residue, 
    following the order of the residues in the PDB (not resi)."""

    # Set color palette
    if not palette:
        palette = palette_dict['rainbow']
    else:
        palette = palette_dict[palette]

    for obj in cmd.get_object_list(f'{selection}'):
        stored.ranked_atoms = {}
        cmd.iterate(obj, "stored.ranked_atoms[rank] = (f'/{model}//{chain}/{resi}', type)")

        # Make a list of residue selections, in the same order as in the PDB
        residues = []
        for k,v in sorted(stored.ranked_atoms.items()):
            if v[1] == 'HETATM':
                continue
            if v[0] not in residues:
                residues.append(v[0])

        # Make color list
        prefix = palette[0]
        min = int(palette[2])
        max = int(palette[3])
        codes = np.linspace(min, max, len(residues), dtype='int')
        colors = ['{}{}'.format(prefix, str(i).zfill(3)) for i in codes]

        # Color residues
        for c, sele in zip(colors, residues):
            cmd.color(c, sele)  
