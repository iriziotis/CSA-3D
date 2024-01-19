from .Superimposer import Superimposer
from .residue_definitions import AA_3TO1, RESIDUE_DEFINITIONS, EQUIVALENT_ATOMS, TEMPLATE_FUZZY_ATOMS
from rmsd import reorder_hungarian


def fit_analogues(p, q, mapping, weighted=False, cycles=1, cutoff=999, scaling_factor=None, reorder=True,
                  transform=False, allow_symmetrics=True, ca=False):
    """Fits two sites on the specified residues (mapping param)."""
    p_atoms, q_atoms, q_coords, q_coords = [],[],[],[]
    for p_res, q_res in mapping.items():
        p_res, q_res = p.residues[p_res], q.residues[q_res]
        _p_atoms, _p_coords = p_res.get_func_atoms()
        _q_atoms, _q_coords = q_res.get_func_atoms()
        p_atoms.extend(_p_atoms)
        p_coords.extend(_p_coords)
        q_atoms.extend(_q_atoms)
        q_coords.extend(_q_coords)




        
e1_id, e2_id = 774, 203
entries_path = '/homes/riziotis/work/research/phd/data/csa3d/entries/entries/'
e1 = pickle.load(open(os.path.join(entries_path, f'csa3d_{str(e1_id).zfill(4)}.ent'), 'rb'))
e2 = pickle.load(open(os.path.join(entries_path, f'csa3d_{str(e2_id).zfill(4)}.ent'), 'rb'))
    
