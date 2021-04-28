import os
import numpy as np
import parity_core
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from .PdbResidue import Het
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
        x,y,z = [],[],[]
        try:
            for res in self.parent_residues:
                # Define boundaries based on functional atoms of residue
                for atom in res.get_atoms():
                    if '{}.{}'.format(res.functional_resname, atom.get_id().upper()) not in RESIDUE_DEFINITIONS:
                        continue
                    x.append(atom.get_coord()[0])
                    y.append(atom.get_coord()[1])
                    z.append(atom.get_coord()[2])
        except (IndexError, KeyError):
            print('Warning: Error occurred while searching for ligands')
            return
        min_x = np.min(x) - headroom
        max_x = np.max(x) + headroom
        min_y = np.min(y) - headroom
        max_y = np.max(y) + headroom
        min_z = np.min(z) - headroom
        max_z = np.max(z) + headroom
        return {'min_x': min_x, 'max_x': max_x,
                'min_y': min_y, 'max_y': max_y,
                'min_z': min_z, 'max_z': max_z}

