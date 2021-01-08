# !/usr/bin/env python3

import rmsd
import numpy as np
from numba import jit, prange
from copy import deepcopy
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.MMCIFParser import FastMMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from .residue_definitions import AA_3TO1, RESIDUE_DEFINITIONS, EQUIVALENT_ATOMS
from .config import COMPOUND_SIMILARITIES, PDB_EC
from .PdbResidue import PdbResidue, Het


class PdbSite:
    """M-CSA PDB catalytic site. A list of PDB objects"""

    def __init__(self):
        self.residues = []
        self.residues_dict = {}
        self.mapped_unisites = []
        self.reference_site = None
        self.structure = None
        self.structure_hets = []
        self.nearby_hets = []
        self.annotations = None

    def __str__(self):
        """Print as pseudo-sequence in one-letter code"""
        return ''.join([AA_3TO1[res.resname] for res in self.residues])

    def __len__(self):
        """Return size of site (residue count)"""
        return self.size

    def __iter__(self):
        """Iterate over residues"""
        yield from self.residues

    def __eq__(self, other):
        """Check if sites contain the same residues"""
        if len(self) == len(other):
            for res in other:
                if res.id not in self.residues_dict:
                    return
            return True
        return

    def __contains__(self, residue):
        """Check if residue of specific id is there"""
        return residue in self.residues

    def __getitem__(self, _id):
        """Return the child with given id."""
        return self.residues_dict[_id]

    def get_residues(self):
        yield from self.residues

    def get_gaps(self):
        """Returns an index of the gap positions (non-aligned residues)"""
        gaps = []
        for i, res in enumerate(self.residues):
            if res.is_gap:
                gaps.append(i)
        return gaps

    def add(self, residue):
        """Add PDB residue to list"""
        if type(residue) == PdbResidue:
            self.residues.append(residue)
            self.residues_dict[residue.id] = residue
            # Add residue structure to site structure
            if residue.structure:
                # Initialise structure
                if self.structure is None:
                    self.structure = Structure(self.id)
                    self.structure.add(Model(0))
                chain_id = residue.structure.get_parent().get_id()
                if chain_id not in self.structure[0]:
                    self.structure[0].add(Chain(chain_id))
                self.structure[0][chain_id].add(residue.structure)
        else:
            print('Attempted to add non-PdbResidue object in PdbSite')
            return
        return True

    def find_ligands(self, parent_structure, headroom=1):
        if type(parent_structure) != Structure:
            return False, False
        all_hets = []
        nearby_hets = []
        residue_list = []
        for residue in self.structure.get_residues():
            chain = residue.get_parent().get_id()
            resid = residue.get_id()[1]
            try:
                residue_list.append(parent_structure[0][chain][resid])
            except KeyError:
                try:
                    for res in parent_structure[0][chain]:
                        if res.get_id()[1] == resid:
                            residue_list.append(res)
                except:
                    print('Warning while searching for ligands:\n'
                          'Could not find residue {} {} from {}'.format(chain, resid, self.pdb_id))
                    continue
        if len(residue_list) == 0:
            return False, False

        # Search for ligands in a box around catalytic residues
        box = Box(residue_list, headroom)
        for residue in parent_structure[0].get_residues():
            residue_id = residue.get_id()
            hetfield = residue_id[0]
            # Capture all HETs in the structure
            if hetfield[0] == 'H':
                het = Het(mcsa_id=self.mcsa_id, pdb_id=self.pdb_id, resname=residue.get_resname(), 
                          resid=residue.get_id()[1], chain=residue.get_parent().get_id())
                het.structure = residue
                all_hets.append(het)
                # Capture HETs in the box
                if het.structure in box:
                    residue.parity_score = Box.similarity_with_cognate(het.structure)
                    residue.centrality = box.mean_distance_from_residues(het.structure)
                    nearby_hets.append(het)
        self.structure_hets = all_hets
        self.nearby_hets = nearby_hets

    def write_pdb(self, write_hets=False, outdir=None, outfile=None):
        """Writes site coordinates in PDB format"""
        if not outdir:
            outdir = '.'
        if not outfile:
            conservation = 'c'
            if not self.is_conserved:
                conservation = 'm'
            elif self.is_conservative_mutation:
                conservation = 'cm'
            outfile = '{}/mcsa_{}.{}.{}.{}.pdb'.format(outdir.strip('/'),
                                                       str(self.mcsa_id).zfill(4), self.id,
                                                       'reference' if self.is_reference else 'cat_site',
                                                       conservation)
        with open(outfile, 'w') as o:
            residues = self.residues.copy()
            if write_hets:
                residues += self.nearby_hets
            if self.ec:
                print('REMARK', self.ec, file=o)
            if self.annotations:
                for k,v in self.annotations.items():
                    print('REMARK', k.upper(), v, file=o)
            for res in residues:
                if res.structure is not None:
                    for atom in res.structure:
                        pdb_line = '{:6}{:5d} {:<4}{}{:>3}{:>2}{:>4}{:>12.3f}' \
                                   '{:>8.3f}{:>8.3f} {:6}'.format(
                            'ATOM' if atom.get_parent().get_id()[0] == ' ' else 'HETATM',
                            atom.get_serial_number() if atom.get_serial_number() else 0,
                            atom.name if len(atom.name) == 4 else ' {}'.format(atom.name),
                            atom.get_altloc(),
                            atom.get_parent().get_resname(),
                            atom.get_parent().get_parent().get_id(),
                            atom.get_parent().get_id()[1],
                            atom.get_coord()[0],
                            atom.get_coord()[1],
                            atom.get_coord()[2],
                            atom.get_occupancy() if atom.get_occupancy() else '')
                        print(pdb_line, file=o)

    def get_atom_strings_and_coords(self, functional_atoms=True, ca_only=False, allow_symmetrics=True):
        """Gets atoms and coordinates for superposition and atom reordering
        calculations

        Args:
            functional_atoms: #TODO

            ca_only: #TODO

            allow_symmetrics: If True, equivalent residues and atoms
                              get the same id string, according to the
                              definitions in residue_definitions.py
                              (EQUIVALENT_ATOMS)
        Returns:
            atoms: A NumPy array of atom identifier strings of type
                   'N.RES.AT' where N is the residue serial number
                   in the .pdb file (consistent among all sites),
                   RES is the residue name and AT is the atom name
            coords: A NumPy array of the atomic coordinates
        """
        atoms = []
        coords = []
        for i, res in enumerate(self.residues):
            # Temporary
            if not res.structure:
                break
            for atom in res.structure:
                resname = res.resname.upper()
                if allow_symmetrics:
                    if 'main' in res.funcloc:
                        resname = 'ANY'
                atmid = '{}.{}'.format(resname, atom.name)
                if atmid in RESIDUE_DEFINITIONS:
                    if allow_symmetrics:
                        if atmid in EQUIVALENT_ATOMS:
                            atmid = EQUIVALENT_ATOMS[atmid]
                    atoms.append('{}.{}'.format(i, atmid))
                    coords.append(atom.get_coord())
        atoms = np.array(atoms)
        coords = np.stack(coords, axis=0)

        return atoms, coords

    def fit(self, other, cycles=10, transform=False, mutate=True, reorder=True, allow_symmetrics=True, get_index=False):
        """Calculates RMSD using the Kabsch algorithm from the rmsd module.
        Can also find the optimal atom alignment using the Hungarian algorithm.
        See https://github.com/charnley/rmsd for more"""
        p = deepcopy(self)
        q = deepcopy(other)

        # TODO check cases where site is fit to reference site.fit(site.reference_structure)
        # I get a stack bug, mcsa 1

        # In case gaps are present, remove the corresponding residues from both structures
        for gap in sorted(q.get_gaps(), reverse=True):
            del p.residues[gap]
            del q.residues[gap]
        # Get atom identifier strings and coords as numpy arrays
        p_atoms, p_coords = p.get_atom_strings_and_coords(allow_symmetrics=allow_symmetrics)
        q_atoms, q_coords = q.get_atom_strings_and_coords(allow_symmetrics=allow_symmetrics)
        q_review = []
        if len(p_atoms) != len(q_atoms):
            raise Exception('Atom number mismatch in sites {} and {}'.format(self.id, other.id))
        # Initial crude superposition
        rot, tran, rms = self._super(p_coords, q_coords, cycles=1)
        q_temp = PdbSite._transform(q_coords, rot, tran)
        # In case of non-conservative mutations, make a pseudo-mutation to
        # facilitate superposition
        if mutate:
            for i, (p_atom, q_atom) in enumerate(zip(p_atoms, q_atoms)):
                if p_atom != q_atom:
                    q_atoms[i] = p_atoms[i]
        # Reorder atoms using the Hungarian algorithm from rmsd package
        if reorder:
            q_review = rmsd.reorder_hungarian(p_atoms, q_atoms, p_coords, q_temp)
            q_coords = q_coords[q_review]
            q_atoms = q_atoms[q_review]
        # Iterative superposition. Get rotation matrix, translation vector and RMSD
        rot, tran, rms = self._super(p_coords, q_coords, cycles, cutoff=6)
        if transform:
            other.structure.transform(rot, tran)
            for het in other.nearby_hets:
                het.structure.transform(rot, tran)
        if get_index:
            return rot, tran, rms, q_review.tolist()
        return rot, tran, rms

    def _super(self, p_coords, q_coords, cycles=10, cutoff=6):
        """Iterative superposition"""
        # TODO Complete docstring
        result = None, None, None
        min_rms = 999
        # Initialize Biopython SVDSuperimposer
        sup = SVDSuperimposer()
        for i in range(cycles):
            sup.set(p_coords, q_coords)
            sup.run()
            rms = sup.get_rms()
            rot, tran = sup.get_rotran()
            if rms < min_rms:
                result = (rot, tran, rms)
            # Transform coordinates
            q_trans = np.dot(q_coords, rot) + tran
            # Find outliers
            diff = np.linalg.norm(p_coords-q_trans, axis=1)
            to_keep = np.where(diff < cutoff)
            # Reject outliers
            p_coords = p_coords[to_keep]
            q_coords = q_coords[to_keep]
        return result

    def _transform(coords, rot, tran):
        """Rotates and translates a set of coordinates (NxD NumPy array)"""
        return np.dot(coords, rot) + tran

    def _map_reference_residues(self):
        """Puts each residue in the site in the correct order, according
        to the reference site, using the individual residue mapping to a
        reference residue. Wherever a mapping cannot be found, an empty
        residue is assigned to that position"""
        if self.reference_site is None:
            return
        for reference_residue in self.reference_site:
            found = False
            for res in self:
                if reference_residue == res.reference_residue:
                    found = True
            if not found:
                gap = PdbResidue()
                gap.reference_residue = reference_residue
                self.add(gap)
        self._reorder()
        return

    def _reorder(self):
        """Residue reordering routine for _map_reference_residues"""
        if self.reference_site is None:
            return
        reorder = []
        for i, reference_residue in enumerate(self.reference_site):
            for j, res in enumerate(self):
                if i == j and reference_residue == res.reference_residue:
                    reorder.append(i)
                elif i != j and reference_residue == res.reference_residue:
                    reorder.append(j)
        self.residues = [self.residues[i] for i in reorder]
        return

    @property
    def mcsa_id(self):
        """Get M-CSA ID of catalytic residues."""
        for res in self.residues:
            if res.mcsa_id:
                return res.mcsa_id
        return

    @property
    def pdb_id(self):
        """Get PDB ID of catalytic residues. Not a unique site ID"""
        for res in self.residues:
            if res.pdb_id:
                return res.pdb_id
        return

    @property
    def ec(self):
        """Get EC number from SIFTS"""
        for res in self.residues:
            if res.chain:
                try:
                    return PDB_EC[(self.pdb_id, res.chain[0])]
                except KeyError:
                    return
        return

    @property
    def id(self):
        """Unique ID of the active site. Consists of PDB ID and a string
        of chain IDs of all residues"""
        return '{}_{}'.format(self.pdb_id, '-'.join(res.chain for res in self.residues))

    @property
    def size(self):
        """Get site size in residue count"""
        return len(self.residues)

    @property
    def is_reference(self):
        """Check if site is reference site"""
        if self.size > 0:
            return self.residues[0].is_reference
        return False

    @property
    def is_conserved(self):
        """Check if all residues are conserved by comparing to the reference"""
        return str(self) == str(self.reference_site) or self.is_reference

    @property
    def is_conservative_mutation(self):
        result = False
        for res in self.residues:
            if not res.is_conserved:
                result = False
            if res.is_conservative_mutation:
                result = True
        return result

    @classmethod
    def from_list(cls, res_list, cif_path=None, annotate=True):
        """Construct PdbSite object directly from residue list"""
        # If we have duplicate residues, each one with a different function
        # location annotation, keep only the one with the side chain
        seen = set()
        to_del = []
        for i, res in enumerate(res_list):
            if res.id in seen:
                to_del.append(i)
            seen.add(res.id)
        for i in (sorted(to_del, reverse=True)):
            if res_list[i].funcloc == '':
                del res_list[i]
        if cif_path:
            parser = FastMMCIFParser(QUIET=True)
            structure = parser.get_structure('', cif_path) 
        site = cls()
        for res in res_list:
            if structure:
                res.add_structure(structure)
            site.add(res)
        if annotate and structure:
            site.annotations = PdbSite.get_annotations(cif_path)
            site.find_ligands(structure)
        return site

    @classmethod
    def build_reference(cls, cat_residues, cif_path=None, annotate=True):
        """Builds reference active site from a list of PDB catalytic residues.
        Assumes that the list only contains one active site, so use it only
        if it is a list of manually annotated catalytic residues"""
        return PdbSite.from_list(cat_residues, cif_path, annotate)

    @classmethod
    def build(cls, seed, reslist, reference_site):
        """Builds active site from a list of catalytic residues that may form
        multiple active sites (e.g. all residues annotated as catalytic in a
        PDB structure). Using a residue as seed, it starts building an active site
        by checking the euclidean distances of all residues that have the same resid
        and name. In the end, it maps the site to the reference defined in the args"""
        site = cls()
        if seed.structure is None:
            return
        for res in reslist:
            candidate = seed.get_nearest_equivalent(res, reslist)
            if candidate is None:
                continue
            if candidate not in site:
                site.add(candidate)
        site.reference_site = reference_site
        site._map_reference_residues()
        return site

    @classmethod
    def build_all(cls, reslist, reference_site, cif_path, annotate=True, redundancy_cutoff=None):
        """Builds all sites in using as input a list of catalytic residues.
        Returns a list of PdbSite objects"""
        sites = []
        # Map structure objects in every residue
        parser = FastMMCIFParser(QUIET=True)
        structure = parser.get_structure('', cif_path) 
        if annotate:
            annotations = PdbSite.get_annotations(cif_path)
        # We want all assembly chains
        new_reslist = []
        for res in reslist:
            for chain in structure[0]:
                if res.chain != chain.get_id()[0]:
                    continue
                try:
                    res_structure = chain[res.auth_resid]
                    new_res = deepcopy(res)
                    new_res.chain = chain.get_id()
                    new_res.structure = res_structure
                    new_reslist.append(new_res)
                except KeyError:
                    continue
        reslist = new_reslist
        # Set a reference residue to make seeds
        ref = None
        for res in reslist:
            if res.auth_resid is not None and res.structure is not None:
                ref = res
                break
        if ref is None:
            return sites
        # Get all equivalents of ref residue and make them seeds
        seeds = []
        for res in reslist:
            if res.is_equivalent(ref):
                if res.structure is None or res in seeds:
                    continue
                seeds.append(res)
        # Build a site from each seed
        for seed in seeds:
            site = cls.build(seed, reslist, reference_site)
            if redundancy_cutoff:
                if len(sites) > 0:
                    _, _, rms = site.fit(sites[-1])
                    if rms < redundancy_cutoff:
                        continue
            if annotate and structure:
                site.find_ligands(structure)
                site.annotations = annotations
            sites.append(site)
        return sites

    @staticmethod
    def get_annotations(cif_path):
        """Parsed an MMCIF file and gets all necessary info for a site.
        Returns a dictionary of annotations"""
        try:
            cif = MMCIF2Dict(cif_path) 
        except:
            return
        annotations = dict()

        annotations['title'] = cif['_struct.title'][0]
        annotations['enzyme'] = cif['_struct.pdbx_descriptor'][0]
        annotations['assembly_id'] = cif['_entity_poly.assembly_id'][0].split('-')[-1]
        annotations['exptl_method'] = cif['_exptl.method']
        annotations['exptl_method'] = annotations['exptl_method'][0]
        if 'nmr' in annotations['exptl_method'].lower():
            annotations['resolution'] = ''
        elif 'microscopy' in annotations['exptl_method'].lower():
            annotations['resolution'] = cif['_em_3d_reconstruction.resolution'][0]
        else:
            try:
                annotations['resolution'] = cif['_refine.ls_d_res_high'][0]
            except:
                annotations['resolution'] = '?'
        try:
            annotations['organism_name'] = cif['_entity_src_nat.pdbx_organism_scientific'][0]
        except KeyError:
            try:
                annotations['organism_name'] = cif['_entity_src_gen.pdbx_gene_src_scientific_name'][0]
            except KeyError:
                annotations['organism_name'] = '?'
        try:
            annotations['organism_id'] = cif['_entity_src_nat.pdbx_ncbi_taxonomy_id'][0]
        except KeyError:
            try:
                annotations['organism_id'] = cif['_entity_src_gen.pdbx_gene_ncbi_taxonomy_id'][0]
            except KeyError:
                annotations['organism_id'] = '?'
        try:
            annotations['host_name'] =  cif['_entity_src_gen.pdbx_host_org_scientific_name'][0]
        except KeyError:
            annotations['host_name'] = '?'[0]
        try:
            annotations['host_id'] = cif['_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id'][0]
        except:
            annotations['host_id'] = '?'

        return annotations


class Box:

    def __init__(self, residue_list, headroom=1):
        """Define a box using the provided residues coordinates adding
        some safe distance around (headroom)"""
        self.residue_list = residue_list
        self.headroom = float(headroom)
        self.points = self._get_boundaries(self.residue_list, self.headroom)
        if not self.points:
            return

    def __contains__(self, residue):
        """Checks if the provided residue is in the box"""
        for atom in residue.get_atoms():
            x = atom.get_coord()[0]
            y = atom.get_coord()[1]
            z = atom.get_coord()[2]
            if self.points['min_x'] <= x <= self.points['max_x'] and \
                    self.points['min_y'] <= y <= self.points['max_y'] and \
                    self.points['min_z'] <= z <= self.points['max_z']:
                return True
        return False

    def _get_boundaries(self, residue_list, headroom):
        """Parse residues (BioPython objects) and get the coordinates
        to make the box"""
        self.coords_x = []
        self.coords_y = []
        self.coords_z = []
        try:
            for residue in self.residue_list:
                for atom in residue:
                    self.coords_x.append(atom.get_coord()[0])
                    self.coords_y.append(atom.get_coord()[1])
                    self.coords_z.append(atom.get_coord()[2])
        except KeyError:
            print('Warning: Error occured while trying to parse structure')
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

    def mean_distance_from_residues(self, component):
        """Calculates the mean distance of an arbitrary residue (protein or HET)
        from the residues that define the box"""
        dist_sum = 0
        nofres = len(self.residue_list)
        for residue in self.residue_list:
            dist_sum += Box.min_distance(residue, component)
        return round(dist_sum / nofres, 3)

    @staticmethod
    def similarity_with_cognate(component):
        """Checks the similarity score of the given compound with the cognate
        ligand of the given pdb, using the PARITY-derived data
        (COMPOUND_SIMILARITIES.json) from Jon"""
        try:
            pdb_id = component.get_full_id()[0]
            hetcode = component.get_resname()
        except:
            return None
        r_key = (pdb_id, hetcode, 'r')
        p_key = (pdb_id, hetcode, 'p')
        if r_key in COMPOUND_SIMILARITIES:
            return float(COMPOUND_SIMILARITIES[r_key])
        elif p_key in COMPOUND_SIMILARITIES:
            return float(COMPOUND_SIMILARITIES[p_key])
        else:
            return None

    @staticmethod
    def min_distance(residue_i, residue_j):
        """Calculates the minimum distance between this residue and residue in argument"""
        distances = []
        for atom_i in residue_i:
            for atom_j in residue_j:
                distances.append(atom_i - atom_j)
        return (min(distances))

