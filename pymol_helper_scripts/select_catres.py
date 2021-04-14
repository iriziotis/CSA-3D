from pymol import cmd, stored
from collections import defaultdict

@cmd.extend
def select_catres(selection):
    """ Makes individual selections, each one containing all residues of the same
    index (in the PDB) from all objects."""

    groups = defaultdict(list)
    for obj in cmd.get_object_list(f'{selection}'):
        stored.ranked_atoms = {}
        cmd.iterate(obj, "stored.ranked_atoms[rank] = f'/{model}//{chain}/{resi}'")

        # Make a list of residue selections, in the same order as in the PDB
        residues = []
        for k,v in sorted(stored.ranked_atoms.items()):
            if v not in residues:
                residues.append(v)

        for i, sele in enumerate(residues):
            groups[i].append(sele)

    # Now make group selections
    for i,reslist in groups.items():
        selection_name = f'catres_{i}'
        cmd.select('temp', 'none')
        for sele in reslist:
            cmd.select('cur', sele)
            cmd.select('temp', '(temp) or (cur)')
        cmd.set_name('temp', selection_name)
    cmd.delete('cur')
