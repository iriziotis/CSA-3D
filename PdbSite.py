import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.MMCIFParser import MMCIFParser, FastMMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rmsd import reorder_hungarian
from .residue_definitions import AA_3TO1, RESIDUE_DEFINITIONS, EQUIVALENT_ATOMS
from .config import COMPOUND_SIMILARITIES, PDB2EC
from .PdbResidue import PdbResidue, Het


class PdbSite:
    """M-CSA PDB catalytic site. Contains lists of PdbResidues and mapped UniProt
    catalytic sites (UniSite objects), a 3D structure (Biopython Structure) build
    from individual PdbResidue structures (Biopython Residue), a parent structure
    (Biopython Structure), all and close to the site hetero components (possible
    ligands according to their chemical similarity to the cognate ligand and their
    centrality in the active site) as Het objects (containing a Biopython Residue
    structure), as well as a dictionary of annotations extracted from the parent
    mmCIF assembly structure file and SIFTS"""

    def __init__(self):
        self.residues = []
        self.residues_dict = {}
        self.mapped_unisites = []
        self.reference_site = None
        self.parent_structure = None
        self.structure = None
        self.structure_hets = []
        self.nearby_hets = []
        self.annotations = None
        # TODO add UniProt IDs from SIFTS -- regardless if we don't have the corresponding sites in M-CSA
        # TODO set reference site of reference as itself

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
        """Check if sites contain the same residues (same IDs)"""
        if len(self) == len(other):
            for res in other:
                if res.id not in self.residues_dict:
                    return False
            return True
        return False

    def __contains__(self, residue):
        """Check if residue is there"""
        return residue in self.residues

    def __getitem__(self, _id):
        """Return the residue with given ID."""
        return self.residues_dict[_id]

    # Alternative constructors

    @classmethod
    def from_list(cls, res_list, cif_path, annotate=True):
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
        site = cls()
        if annotate:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure('', cif_path)
            annotations = PdbSite._get_annotations(parser._mmcif_dict)
        else:
            parser = FastMMCIFParser(QUIET=True)
            structure = parser.get_structure('', cif_path)
        for res in res_list:
            if structure:
                res.add_structure(structure)
            site.add(res)
        if annotate:
            site.annotations = annotations
            site.find_ligands()
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
        # Map structure objects in every residue
        if annotate:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure('', cif_path) 
            annotations = PdbSite._get_annotations(parser._mmcif_dict)
        else:
            parser = FastMMCIFParser(QUIET=True)
            structure = parser.get_structure('', cif_path) 
        # We want all equivalent residues from identical assembly chains
        reslist = PdbSite._get_assembly_residues(reslist, structure)
        sites = []
        # Set a reference residue to make seeds
        seeds = PdbSite._get_seeds(reslist)
        # Build a site from each seed
        for seed in seeds:
            site = cls.build(seed, reslist, reference_site)
            if site.has_missing_functional_atoms:
                continue
            if redundancy_cutoff:
                if len(sites) > 0:
                    _, _, _, rms_all = site.fit(sites[-1])
                    if rms_all < redundancy_cutoff:
                        continue
            if annotate and structure:
                site.parent_structure = structure
                site.annotations = annotations
                site.find_ligands()
            sites.append(site)
        return sites

    # Properties

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
                    return PDB2EC[(self.pdb_id, res.chain[0])]
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
    def is_conservative_mutation(self, ignore_funcloc_main=True):
        result = False
        for res in self.residues:
            if ignore_funcloc_main:
                if 'main' in res.funcloc:
                    result = True
                    continue
            if not res.is_conserved and not res.is_conservative_mutation:
                return False
            if res.is_conservative_mutation:
                result = True
        return result

    @property
    def has_missing_functional_atoms(self):
        gaps = set(self.get_gaps())
        self_atoms, _ = self._get_atom_strings_and_coords(omit=gaps)
        ref_atoms, _ = self.reference_site._get_atom_strings_and_coords(omit=gaps)
        return len(self_atoms) != len(ref_atoms)

    def add(self, residue):
        """Add PdbResidue object to site (in the residues list and dict)"""
        if type(residue) == PdbResidue:
            self.residues.append(residue)
            self.residues_dict[residue.id] = residue
            # Add residue structure to site structure
            if residue.structure:
                # Initialize structure if empty
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

    def get_distances(self):
        """Calculates all intra-site residue distances and returns a
        numpy array"""
        dists = []
        seen = set()
        for p in self.residues:
            for q in self.residues:
                if p==q or (q.id, p.id) in seen or p.is_gap or q.is_gap:
                    continue
                dists.append(p-q)
                seen.add((p.id, q.id))
        return np.array(dists)

    def get_residues(self):
        """To iterate over catalytic residues"""
        yield from self.residues

    def get_gaps(self):
        """Returns an index of the gap positions (non-aligned residues)"""
        gaps = []
        for i, res in enumerate(self.residues):
            if res.is_gap:
                gaps.append(i)
        return gaps

    def find_ligands(self, headroom=1):
        """
        Searches the parent structure for hetero components close to the
        catalytic residues, by defining a box area around them plus some
        extra distance. Populates the nearby_hets list with Het objects

        Args:
            headroom: the extra search space (in Å) around the catalytic residues
        """
        if type(self.parent_structure) != Structure:
            return False, False
        all_hets = []
        nearby_hets = []
        residue_list = []
        for residue in self.structure.get_residues():
            chain = residue.get_parent().get_id()
            resid = residue.get_id()[1]
            try:
                residue_list.append(self.parent_structure[0][chain][resid])
            except (IndexError, KeyError):
                try:
                    for res in self.parent_structure[0][chain]:
                        if res.get_id()[1] == resid:
                            residue_list.append(res)
                except (IndexError, KeyError):
                    continue
        if len(residue_list) == 0:
            return False, False
        # Search for ligands in a box around catalytic residues
        box = Box(residue_list, headroom)
        for residue in self.parent_structure[0].get_residues():
            residue_id = residue.get_id()
            hetfield = residue_id[0]
            # Capture all HETs in the structure
            if hetfield[0] == 'H' and 'HOH' not in hetfield:
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

    # TODO write ligands as REMARK entries
    def write_pdb(self, write_hets=False, outdir=None, outfile=None):
        """
        Writes site coordinates in PDB format
        Args:
            write_hets: Include coordinates of nearby hets.
            outdir: Directory to save the .pdb file
            outfile: If unspecified, name is formatted to include
                     info on M-CSA ID, chain of each catalytic residue,
                     annotation if the site is a reference site and
                     an annotation about the conservation, relatively
                     to the reference (c: conserved, m: mutated, cm: has
                     only conservative mutations)
        """
        if not outdir:
            outdir = '.'
        if not outfile:
            conservation = 'c'
            if not self.is_conserved:
                conservation = 'm'
            elif self.is_conservative_mutation:
                conservation = 'cm'
            outfile = '{}/mcsa_{}.{}.{}.{}.pdb'.format(outdir.strip('/'), str(self.mcsa_id).zfill(4), self.id,
                      'reference' if self.is_reference else 'cat_site', conservation)
        with open(outfile, 'w') as o:
            residues = self.residues.copy()
            if write_hets:
                residues += self.nearby_hets
            if self.ec:
                print('REMARK', self.ec, file=o)
            if self.annotations:
                for k, v in self.annotations.items():
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
            print('END', file=o)

    def fit(self, other, cycles=10, transform=False, mutate=True, reorder=True, allow_symmetrics=True):
        """Iteratively fits two catalytic sites (self: fixed site, other: mobile site)
        using the Kabsch algorithm from the rmsd module (https://github.com/charnley/rmsd).
        Can also find the optimal atom alignment in each residue, considering
        symmetrical atoms and functionally similar residues, using the
        Hungarian algorithm.

        Args:
            other: mobile active site to fit
            cycles: Number of fitting iterations to exclude outlying atoms
            transform: Also transforms the mobile site's coordinates
            mutate: If the two active sites do not have the same residues,
                    make pseudo-mutations to the mobile site to facilitate
                    atom correspondence
            reorder: Find the optimal atom correspondence (within a residue)
                     between the two sites, taking into account conservative
                     mutations and symmetrical atoms (optional). See and
                     definitions in residue_definitions.py module.
            allow_symmetrics: Allows flipping of side chains if atoms are
                              equivalent or symmetrical
            outlier_rms: TODO compute rms over all functional atoms (including the ones
                         that where excluded during the outlier rejection
            get_index: To return the new order of atoms after applying the
                       Hungarian algorithm as an index list.

        Returns: rot, tran, rms, rms_all
            rot: Rotation matrix to transform mobile site into the fixed site
            tran: Translation vector to transform mobile site into the fixed site
            rms: RMSD after fitting, excluding outliers
            rms_all: RMSD over all atoms, including outliers

        Raises:
            Exception: If number of functions atoms in the two sites is not the same (e.g.
                       if there are missing atoms from the parent structure)
        """
        # In case gaps are present, exclude those positions
        gaps = set(self.get_gaps() + other.get_gaps())
        # Get atom identifier strings and coords as numpy arrays
        p_atoms, p_coords = self._get_atom_strings_and_coords(allow_symmetrics, omit=gaps)
        q_atoms, q_coords = other._get_atom_strings_and_coords(allow_symmetrics, omit=gaps)
        if len(p_atoms) != len(q_atoms):
            raise Exception('Atom number mismatch in sites {} and {}'.format(self.id, other.id))
        # Initial crude superposition
        rot, tran, rms, _ = PdbSite._super(p_coords, q_coords, cycles=1)
        q_trans = PdbSite._transform(q_coords, rot, tran)
        # In case of non-conservative mutations, make pseudo-mutations to facilitate superposition
        if mutate:
            for i, (p_atom, q_atom) in enumerate(zip(p_atoms, q_atoms)):
                if p_atom != q_atom:
                    q_atoms[i] = p_atoms[i]
        # Reorder atoms using the Hungarian algorithm from rmsd package
        if reorder:
            q_review = reorder_hungarian(p_atoms, q_atoms, p_coords, q_trans)
            q_coords = q_coords[q_review]
        # Iterative superposition. Get rotation matrix, translation vector and RMSD
        rot, tran, rms, rms_all = PdbSite._super(p_coords, q_coords, cycles, cutoff=6)
        if transform:
            other.structure.transform(rot, tran)
            for het in other.nearby_hets:
                het.structure.transform(rot, tran)
        return rot, tran, rms, rms_all

    # Private methods

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

    def _get_atom_strings_and_coords(self, functional_atoms=True, ca_only=False, allow_symmetrics=True, omit=None):
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
            if omit:
                if i in omit:
                    continue
            if not res.structure:
                return
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

    @staticmethod
    def _get_assembly_residues(reslist, parent_structure):
        """
        Makes a new residue list of all equivalent residues found in identical assembly
        chains. Also applies an auth_resid correction  where residues in identical chains 
        might have a different auth_resid (usually of 1xxx or 2xxx for chains A and B 
        respectively). Structures are also mapped in the residues.

        Args:
            reslist: The residue list to be enriched
            parent_structure: BioPython Structure object of the parent structure
        Returns:
            An enriched list of residues with mapped structures.
        """
        new_reslist = []
        for res in reslist:
            for chain in parent_structure[0]:
                if res.chain != chain.get_id()[0]:
                    continue
                try:
                    res_structure = chain[res.auth_resid]
                except KeyError:
                    try: 
                        res_structure = chain[res.corrected_auth_resid]
                    except KeyError:
                        continue
                new_res = res.copy()
                new_res.chain = chain.get_id()
                new_res.structure = res_structure
                new_reslist.append(new_res)
        return new_reslist
    
    @staticmethod
    def _get_seeds(reslist):
        """Finds residues in a list of that can be used as seeds when
        building multiple active sites"""
        seeds = []
        ref = None
        for res in reslist:
            if res.auth_resid is not None and res.structure is not None:
                ref = res
                break
        if ref is None:
            return seeds
        # Get all equivalents of ref residue and make them seeds
        for res in reslist:
            if res.is_equivalent(ref):
                if res.structure is None or res in seeds:
                    continue
                seeds.append(res)
        return seeds

    @staticmethod
    def _super(p_coords, q_coords, cycles=10, cutoff=6):
        """
        Iterative superposition of two coordinate sets (NumPy arrays)

        Args:
            p_coords: Fixed coord set
            q_coords: Mobiles coord set
            cycles: Number of outlier rejection iterations
            cutoff: Pairwise atom distance threshold to reject outliers

        Returns: rot, tran, rms, rms_all
            rot: Rotation matrix
            tran: Translation vector
            rms: RMSD after fitting (not including outliers)
            rms_all: RMSD over all atoms, including outliers
        """
        min_rms = 999
        p_all = np.array(p_coords, copy=True)
        q_all = np.array(q_coords, copy=True)
        # Initialize Biopython SVDSuperimposer
        sup = SVDSuperimposer()
        for i in range(cycles):
            sup.set(p_coords, q_coords)
            sup.run()
            rms = sup.get_rms()
            rot, tran = sup.get_rotran()
            if rms < min_rms:
                results = (rot, tran, rms)
            # Transform coordinates
            q_trans = np.dot(q_coords, rot) + tran
            # Find outliers
            diff = np.linalg.norm(p_coords-q_trans, axis=1)
            to_keep = np.where(diff < cutoff)
            # Reject outliers
            p_coords = p_coords[to_keep]
            q_coords = q_coords[to_keep]
        # Also compute RMSD over all atoms
        rot, tran, rms = results
        q_all = np.dot(q_all, rot) + tran
        rms_all = PdbSite._rmsd(p_all, q_all)
        return rot, tran, np.round(rms, 3), np.round(rms_all, 3)

    @staticmethod
    def _rmsd(p_coords, q_coords):
        """Calculates rmsd on two coordinate sets (NumPy arrays) WITHOUT
        transformation and minimization"""
        diff = np.square(np.linalg.norm(p_coords-q_coords, axis=1))
        return np.sqrt(np.sum(diff)/diff.size)

    @staticmethod
    def _transform(coords, rot, tran):
        """Rotates and translates a set of coordinates (NxD NumPy array)"""
        return np.dot(coords, rot) + tran

    @staticmethod
    def _get_annotations(cif):
        """
        Parses an MMCIF raw file or dictionary and gets all necessary info for a site.

        Args:
            cif: Either the path of the raw mmCIF file, or a pre-made dictionary from
                 MMCIF2Dict()
        Returns: A smaller dictionary of annotations with less complex key names
        """
        if not isinstance(cif, dict):
            if not isinstance(cif, str) or not cif_path.endswith('cif'):
                return
            cif = MMCIF2Dict(cif_path)
        annotations = dict()
        annotations['title'] = cif['_struct.title'][0]
        annotations['enzyme'] = cif['_struct.pdbx_descriptor'][0]
        annotations['assembly_id'] = cif['_entity_poly.assembly_id'][0].split('-')[-1]
        annotations['exptl_method'] = cif['_exptl.method'][0]
        if 'nmr' in annotations['exptl_method'].lower():
            annotations['resolution'] = ''
        elif 'microscopy' in annotations['exptl_method'].lower():
            annotations['resolution'] = cif['_em_3d_reconstruction.resolution'][0]
        else:
            try:
                annotations['resolution'] = cif['_refine.ls_d_res_high'][0]
            except KeyError:
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
            annotations['host_name'] = cif['_entity_src_gen.pdbx_host_org_scientific_name'][0]
        except KeyError:
            annotations['host_name'] = '?'[0]
        try:
            annotations['host_id'] = cif['_entity_src_gen.pdbx_host_org_ncbi_taxonomy_id'][0]
        except KeyError:
            annotations['host_id'] = '?'

        return annotations


class Box:

    def __init__(self, residue_list, headroom=1):
        """Define a box using the provided residues coordinates adding
        some safe distance around (headroom)"""
        self.residue_list = residue_list
        self.headroom = float(headroom)
        self.points = self._get_boundaries(self.headroom)
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

    def _get_boundaries(self, headroom):
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
        except IndexError:
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
        return min(distances)
