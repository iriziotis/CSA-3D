import os
import numpy as np
import parity_core
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from .PdbResidue import PdbResidue, Het
from .config import WORKING_DIR, REACTION_MOLS_DIR, HET_MOLS_DIR,  COMPOUND_SIMILARITIES, PDB2EC, PDB2UNI, HET2SMILES, EC_REACTION
from .residue_definitions import RESIDUE_DEFINITIONS

class LigandBox:

    def __init__(self, site, headroom=1):
        """Define a box-like area using the functional atoms of the site, 
        plus some safe distance around (headroom)"""
        self.site = site
        self.headroom = float(headroom)
        self.parent_residues = self._get_parent_residues()
        self.points = self._get_boundaries(self.headroom)

    def __contains__(self, entity):
        """Checks if the provided het (or residue, or atom) is in the box"""
        if type(entity) == Het:
            atom_list = entity.structure.get_atoms()
        elif type(entity) == Residue:
            atom_list = entity.get_atoms()
        for atom in atom_list:
            x,y,z = atom.get_coord()
            if self.points['min_x'] <= x <= self.points['max_x'] and \
                    self.points['min_y'] <= y <= self.points['max_y'] and \
                    self.points['min_z'] <= z <= self.points['max_z']:
                return True
        return False

    #def mean_distance_from_residues(self, het):
    #    """Calculates the mean distance of an arbitrary residue (protein or HET)
    #    from the residues that define the box"""
    #    dist_sum = 0
    #    nofres = len(self.parent_residues)
    #    for residue in self.parent_residues:
    #        dist_sum += LigandBox.min_distance(residue, het.structure)
    #    return round(dist_sum / nofres, 3)

    def mean_distance_from_residues(self, het):
        dists = [] 
        for res in self.site:
            dists.append(het.get_distance(res, kind='com'))
        return np.nanmean(np.array(dists))

    def similarity_with_cognate(self, het):
        try:
            pdb_id = het.pdb_id
            hetcode = het.resname.upper()
        except IndexError:
            return None
        if het.flag != 'P':
            try:
                return self.generate_parity(hetcode)
            except Exception as e:
                return None
        return None
        
    #def similarity_with_cognate(self, het):
    #    """Checks the similarity score of the given compound with the cognate
    #    ligand of the given pdb, using the PARITY-derived data"""
    #    try:
    #        pdb_id = het.pdb_id
    #        hetcode = het.resname.upper()
    #    except IndexError:
    #        return None
    #    r_key = (pdb_id, hetcode, 'r')
    #    p_key = (pdb_id, hetcode, 'p')
    #    if r_key in COMPOUND_SIMILARITIES:
    #        return float(COMPOUND_SIMILARITIES[r_key])
    #    elif p_key in COMPOUND_SIMILARITIES:
    #        return float(COMPOUND_SIMILARITIES[p_key])
    #    else:
    #        if het.flag != 'P':
    #            try:
    #                return self.generate_parity(hetcode)
    #            except Exception as e:
    #                return None

    def generate_parity(self, hetcode):
        ec = self.site.ec
        try:
            b_sdf = f'{HET_MOLS_DIR}/{hetcode}.sdf'
            cognate_reactants, cognate_products = EC_REACTION[ec]
        except KeyError:
            return
        max_score = 0.0
        match = None
        component_type = None
        for components in (cognate_reactants, cognate_products):
            for c in components:
                c_sdf = f'{REACTION_MOLS_DIR}/{c}.mol'
                if not os.path.exists(c_sdf) or not os.path.exists(b_sdf):
                    continue
                try:
                    score = parity_core.generate_parity(b_sdf, c_sdf)
                except Exception as e:
                    continue
                if score >= max_score:
                    max_score = score
                    match = c
                    component_type = 'r' if components is cognate_reactants else 'p'

        # Just to compile a dataset, will be removed afterwards
        line = '{},{},{},{},{},{}'.format(self.site.pdb_id, ec, hetcode, match, component_type, np.round(max_score, 3))
        os.makedirs(f'{WORKING_DIR}/parity_all', exist_ok=True)
        outfile = '{}/parity_all/csa3d_{}.parity_all.csv'.format(WORKING_DIR, str(self.site.mcsa_id).zfill(4))
        if os.path.exists(outfile):
            if line in open(outfile).read():
                return max_score
        with open(outfile, 'a') as o:
            print(line, file=o)

        return max_score

    @staticmethod
    def min_distance(residue_i, residue_j):
        """Calculates the minimum distance between this residue and residue in argument"""
        distances = []
        for atom_i in residue_i:
            for atom_j in residue_j:
                distances.append(atom_i - atom_j)
        return min(distances)

    def _get_parent_residues(self):
        """Populates the list of the corresponding residues from parent structure (untransformed)"""
        parent_residues = []
        for residue in self.site:
            try:
                chain = residue.structure.get_parent().get_id()
                resid = residue.structure.get_id()[1]
            except AttributeError:
                continue
            try:
                parent_residue = self.site.parent_structure[0][chain][resid]
            except (IndexError, KeyError):
                try:
                    for res in self.site.parent_structure[0][chain]:
                        if res.get_id()[1] == resid:
                            parent_residue = res
                except (IndexError, KeyError):
                    continue
            resname = residue.resname.upper()
            if residue.has_main_chain_function or not residue.is_standard:
                resname = 'ANY'
            parent_residue.functional_resname = resname
            parent_residues.append(parent_residue)
        return parent_residues

    def _get_boundaries(self, headroom):
        """Parse residues (BioPython objects) and get the coordinates
        to make the box"""
        self.coords_x = []
        self.coords_y = []
        self.coords_z = []
        try:
            for res in self.parent_residues:
                # Define boundaries based on functional atoms of residue
                for atom in res.get_atoms():
                    if '{}.{}'.format(res.functional_resname, atom.get_id().upper()) not in RESIDUE_DEFINITIONS:
                        continue
                    self.coords_x.append(atom.get_coord()[0])
                    self.coords_y.append(atom.get_coord()[1])
                    self.coords_z.append(atom.get_coord()[2])
        except (IndexError, KeyError):
            print('Warning: Error occurred while trying to parse structure')
            return
        min_x = min(self.coords_x) - headroom
        max_x = max(self.coords_x) + headroom
        min_y = min(self.coords_y) - headroom
        max_y = max(self.coords_y) + headroom
        min_z = min(self.coords_z) - headroom
        max_z = max(self.coords_z) + headroom
        return {'min_x': min_x, 'max_x': max_x,
                'min_y': min_y, 'max_y': max_y,
                'min_z': min_z, 'max_z': max_z}
