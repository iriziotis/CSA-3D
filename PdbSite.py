import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.MMCIFParser import MMCIFParser, FastMMCIFParser
from rmsd import reorder_hungarian
from .residue_definitions import AA_3TO1, RESIDUE_DEFINITIONS, EQUIVALENT_ATOMS
from .config import COMPOUND_SIMILARITIES, PDB2EC, PDB2UNI
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
        self.mmcif_dict = dict()

    def __str__(self):
        """Show as pseudo-sequence in one-letter code"""
        return self.sequence

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
        return residue.full_id in self.residues_dict

    def __getitem__(self, full_id):
        """Return the residue with given ID."""
        return self.residues_dict[full_id]

    # Alternative constructors

    @classmethod
    def from_list(cls, reslist, cif_path, annotate=True):
        """Construct PdbSite object directly from residue list"""
        # If we have duplicate residues, each one with a different function
        # location annotation, keep only the one with the side chain
        mmcif_dict = dict()
        # First reduce redundant residues with multiple function locations
        reslist = PdbSite._cleanup_list(reslist)
        site = cls()
        if annotate:
            try:
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure('', cif_path)
                mmcif_dict = parser._mmcif_dict
            except TypeError:
                return
        else:
            parser = FastMMCIFParser(QUIET=True)
            structure = parser.get_structure('', cif_path)
        for res in reslist:
            if structure:
                res.add_structure(structure)
            site.add(res)
        if annotate:
            site.mmcif_dict = mmcif_dict
            site.find_ligands()
        return site

    @classmethod
    def build_reference(cls, reslist, cif_path, annotate=True):
        """Builds reference active site from a list of PDB catalytic residues.
        Assumes that the list only contains one active site, so use it only
        if it is a list of manually annotated catalytic residues"""
        ref = PdbSite.from_list(reslist, cif_path, annotate)
        ref.reference_site = ref
        return ref

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
            candidate = PdbSite._get_nearest_equivalent(res, seed, reslist, site)
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
        sites = []
        mmcif_dict = dict()
        try:
            if annotate:
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure('', cif_path)
                mmcif_dict = parser._mmcif_dict
            else:
                parser = FastMMCIFParser(QUIET=True)
                structure = parser.get_structure('', cif_path)
        except TypeError:
            return sites
        # First reduce redundant residues with multiple function locations
        reslist = PdbSite._cleanup_list(reslist)
        # We want all equivalent residues from identical assembly chains
        reslist = PdbSite._get_assembly_residues(reslist, structure)
        # Get seeds to build active sites
        seeds = PdbSite._get_seeds(reslist)
        # Build a site from each seed
        for seed in seeds:
            site = cls.build(seed, reslist, reference_site)
            if site.has_missing_functional_atoms:
                continue
            # Reduce redundancy
            if len(sites) > 0:
                if site.has_identical_residues(sites[-1]):
                    continue
                if redundancy_cutoff:
                    _, _, _, rms_all = site.fit(sites[-1])
                    if rms_all < redundancy_cutoff:
                        continue
            # Add ligands and annotations
            if annotate and structure:
                site.parent_structure = structure
                site.mmcif_dict = mmcif_dict
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
    def uniprot_id(self):
        """Get UniProt ID of the chain of the first residue"""
        for res in self.residues:
            if res.chain:
                try:
                    return PDB2UNI[(self.pdb_id, res.chain[0])]
                except KeyError:
                    continue
        return

    @property
    def ec(self):
        """Get EC number from SIFTS"""
        for res in self.residues:
            if res.chain:
                try:
                    return PDB2EC[(self.pdb_id, res.chain[0])]
                except KeyError:
                    continue
        return

    @property
    def sequence(self):
        """Show as pseudo-sequence in one-letter code"""
        return ''.join([AA_3TO1[res.resname] for res in self.residues])

    @property
    def title(self):
        """Return title of PDB entry"""
        try:
            return self.mmcif_dict['_struct.title'][0]
        except KeyError:
            return

    @property
    def enzyme(self):
        """Return enzyme name"""
        try:
            return self.mmcif_dict['_struct.pdbx_descriptor'][0]
        except KeyError:
            return

    @property
    def assembly_id(self):
        """Return PDB assembly ID"""
        try:
            return int(self.mmcif_dict['_entity_poly.assembly_id'][0][-1])
        except (TypeError, KeyError):
            return

    @property
    def experimental_method(self):
        """Return structure determination method"""
        try:
            return self.mmcif_dict['_exptl.method'][0]
        except KeyError:
            return

    @property
    def resolution(self):
        """Return resolution in Angstrom"""
        try:
            if 'nmr' in self.experimental_method.lower():
                return
            elif 'microscopy' in self.experimental_method.lower():
                return float(self.mmcif_dict['_em_3d_reconstruction.resolution'][0])
            else:
                return float(self.mmcif_dict['_refine.ls_d_res_high'][0])
        except (TypeError, KeyError, AttributeError):
            return

    @property
    def organism_name(self):
        """Return name of organism of origin"""
        try:
            return self.mmcif_dict['_entity_src_nat.pdbx_organism_scientific'][0]
        except KeyError:
            try:
                return self.mmcif_dict['_entity_src_gen.pdbx_gene_src_scientific_name'][0]
            except KeyError:
                return

    @property
    def organism_id(self):
        """Return id of organism of origin"""
        try:
            return self.mmcif_dict['_entity_src_nat.pdbx_ncbi_taxonomy_id'][0]
        except KeyError:
            try:
                return self.mmcif_dict['_entity_src_gen.pdbx_gene_src_ncbi_taxonomy_id'][0]
            except KeyError:
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
        """Checks if the mutations in the site are conservative. Option to
        ignore residues that function via main chain"""
        result = False
        for res in self.residues:
            if ignore_funcloc_main:
                if res.has_main_chain_function or res.has_double_funcloc:
                    result = True
                    continue
            if not res.is_conserved and not res.is_conservative_mutation:
                return False
            if res.is_conservative_mutation:
                result = True
        return result

    @property
    def has_missing_functional_atoms(self):
        """Checks if there are missing functional atoms from the residue
        structures or site is empty"""
        try:
            gaps = set(self.get_gaps())
            self_atoms, _ = self._get_atom_strings_and_coords(omit=gaps)
            ref_atoms, _ = self.reference_site._get_atom_strings_and_coords(omit=gaps)
            return len(self_atoms) != len(ref_atoms)
        except ValueError:
            return True

    # Methods

    def add(self, residue):
        """Add PdbResidue object to site (in the residues list and dict)"""
        if type(residue) == PdbResidue:
            self.residues.append(residue)
            self.residues_dict[residue.full_id] = residue
            if residue.structure:
                # Initialize structure if empty
                if self.structure is None:
                    self.structure = Structure(self.id)
                    self.structure.add(Model(0))
                chain_id = residue.structure.get_parent().get_id()
                if chain_id not in self.structure[0]:
                    self.structure[0].add(Chain(chain_id))
                # Add residue structure to site structure
                self.structure[0][chain_id].add(residue.structure)
        else:
            print('Attempted to add non-PdbResidue object in PdbSite')
            return False
        return True

    def get_distances(self):
        """Calculates all intra-site residue distances and returns a
        numpy array"""
        dists = []
        seen = set()
        for p in self.residues:
            for q in self.residues:
                if p == q or (q.id, p.id) in seen or p.is_gap or q.is_gap:
                    continue
                dists.append(p - q)
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

    def contains_equivalent(self, res):
        """Checks if the site contains a catalytic residue of the basic info
        (name, resid, auth_resid), and either the same chiral_id or chain"""
        for sres in self:
            if sres.is_equivalent(res, by_chiral_id=True) or \
               sres.is_equivalent(res, by_chiral_id=False, by_chain=True):
                return True
        return False

    def has_identical_residues(self, other):
        """Checks if two sites have the same residues, although their order might be
        different. Used to cleanup redundant symmetrical active sites like
        HIV-protease"""
        for res in other:
            if not self.contains_equivalent(res):
                return False
        return True

    def get_chiral_residues(self):
        """Gets chiral residues from the site if there are any (residues that have
        the same resname, resid, auth_resid but different chains)"""
        identicals = set()
        for i in self:
            for j in self:
                if i==j or i.is_gap or j.is_gap:
                    continue
                if i.is_equivalent(j, by_chiral_id=False, by_chain=False):
                    if (j.chiral_id, i.chiral_id) not in identicals:
                        identicals.add((i.chiral_id, j.chiral_id))
        return identicals

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
                if het in box:
                    het.parity_score = box.similarity_with_cognate(het)
                    het.centrality = box.mean_distance_from_residues(het)
                    nearby_hets.append(het)
        self.structure_hets = all_hets
        self.nearby_hets = nearby_hets

    def write_pdb(self, outdir=None, outfile=None, write_hets=False, func_atoms_only=False):
        """
        Writes site coordinates in PDB format
        Args:
            write_hets: Include coordinates of nearby hets.
            outdir: Directory to save the .pdb file
            outfile: If unspecified, name is formatted to include info on M-CSA ID, 
                     chain of each catalytic residue, annotation if the site is a
                     reference site and an annotation about the conservation, relatively
                     to the reference (c: conserved, m: mutated, cm: has only conservative
                     mutations)
        """
        if not outdir:
            outdir = '.'
        if not outfile:
            conservation = 'm'
            if self.is_conservative_mutation:
                conservation = 'cm'
            if self.is_conserved:
                conservation = 'c'
            if func_atoms_only:
                atms = 'func'
            else:
                atms = 'all'
            outfile = '{}/mcsa_{}.{}.{}.{}.{}.pdb'.format(outdir.strip('/'), str(self.mcsa_id).zfill(4), self.id,
                                                       'reference' if self.is_reference else 'cat_site', conservation, atms)
        with open(outfile, 'w') as o:
            if bool(self.mmcif_dict):
                all_hets = ','.join('{0.resname};{0.resid};{0.parity_score};{0.centrality}'.format(h) for h in self.structure_hets)
                nearby_hets = ','.join('{0.resname};{0.resid};{0.parity_score};{0.centrality}'.format(h) for h in self.nearby_hets)
                remarks = ('REMARK CATALYTIC SITE\n'
                           'REMARK ID {0.id}\n'
                           'REMARK PDB_ID {0.pdb_id}\n'
                           'REMARK ASSEMBLY_ID {0.assembly_id}\n'
                           'REMARK UNIPROT_ID {0.uniprot_id}\n'
                           'REMARK EC {0.ec}\n'
                           'REMARK TITLE {0.title}\n'
                           'REMARK ENZYME {0.enzyme}\n'
                           'REMARK EXPERIMENTAL_METHOD {0.experimental_method}\n'
                           'REMARK RESOLUTION {0.resolution}\n'
                           'REMARK ORGANISM_NAME {0.organism_name}\n'
                           'REMARK ORGANISM_ID {0.organism_id}\n'
                           'REMARK ALL_HETS {1}\n'
                           'REMARK NEARBY_HETS {2}'.format(self, all_hets, nearby_hets))
                print(remarks, file=o)
            residues = self.residues.copy()
            if write_hets:
                residues += self.nearby_hets
            for res in residues:
                if res.structure is not None:
                    for atom in res.structure:
                        if func_atoms_only and type(res) == PdbResidue:
                            resname = res.resname.upper()
                            if res.has_main_chain_function or not res.is_standard:
                                resname = 'ANY'
                            if '{}.{}'.format(resname, atom.get_id().upper()) not in RESIDUE_DEFINITIONS:
                                continue
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

    def fit(self, other, cycles=10, cutoff=6, transform=False, mutate=True, reorder=True, allow_symmetrics=True):
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
        rot, tran, rms, rms_all = PdbSite._super(p_coords, q_coords, cycles, cutoff)
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
                gap = PdbResidue(mcsa_id=reference_residue.mcsa_id, 
                                 pdb_id=reference_residue.pdb_id, 
                                 chiral_id=reference_residue.chiral_id)
                gap.reference_residue = reference_residue
                self.add(gap)
        self._reorder()
        return

    def _get_atom_strings_and_coords(self, allow_symmetrics=True, omit=None):
        """Gets atoms and coordinates for superposition and atom reordering
        calculations

        Args:
            allow_symmetrics: If True, equivalent residues and atoms
                              get the same id string, according to the
                              definitions in residue_definitions.py
                              (EQUIVALENT_ATOMS)
            omit: Residues to exclude
        Returns:
            atoms: A NumPy array of atom identifier strings of type
                   'N.RES.AT' where N is the residue serial number
                   in the .pdb file (consistent among all sites),
                   RES is the residue name and AT is the atom name
            coords: A NumPy array of the atomic coordinates
        """
        atoms = []
        coords = []
        for i, res in enumerate(self):
            if omit:
                if i in omit:
                    continue
            if not res.structure:
                return atoms, coords
            for atom in res.structure:
                resname = res.resname.upper()
                if allow_symmetrics:
                    if res.has_main_chain_function:
                        resname = 'ANY'
                    if not res.is_standard:
                        resname = 'PTM'
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

        # If site contains chiral residues, reorder them by chain
        chiral = self.get_chiral_residues()
        if chiral:
            for pair in chiral:
                p = self.residues[pair[0]]
                q = self.residues[pair[1]]
                if q.chain < p.chain:
                    self.residues[pair[0]], self.residues[pair[1]] = \
                    self.residues[pair[1]], self.residues[pair[0]]
        return

    @staticmethod
    def _cleanup_list(reslist):
        """Finds duplicate residues of different funclocs and makes a single
        one with two funclocs. Returns a new list without redundant residues"""
        new_reslist = []
        seen = set()
        ignore = set()
        for p in reslist:
            for q in reslist:
                if p==q or (q.id, p.id) in seen:
                    continue
                if p.is_equivalent(q, by_chiral_id=False, by_chain=True):
                    if p.funclocs != q.funclocs:
                        new_res = p.copy()
                        new_res.funclocs = [p.funclocs[0], q.funclocs[0]]
                        new_reslist.append(new_res)
                        ignore.add(p.id)
                        ignore.add(q.id)
                seen.add((p.id, q.id))
        for p in reslist:
            if p.id not in ignore and p not in new_reslist:
                new_reslist.append(p)
        return new_reslist

    @staticmethod
    def _get_assembly_residues(reslist, parent_structure):
        """
        Makes a new residue list of all equivalent residues found in identical assembly
        chains. Also applies an auth_resid correction  where residues in identical chains 
        might have a different auth_resid (usually of 1xxx or 2xxx for chains A and B 

        Args:
            reslist: The residue list to be enriched
            parent_structure: BioPython Structure object of the parent structure
        Returns:
            An enriched list of residues with mapped structures.
        """
        new_reslist = []
        for res in reslist:
            res_structure = None
            for chain in parent_structure[0]:
                if res.chain != chain.get_id()[0]:
                    continue
                # If we have a standard residue
                if res.is_standard:
                    try:
                        res_structure = chain[res.auth_resid]
                    except KeyError:
                        try:
                            res_structure = chain[res.corrected_auth_resid]
                        except KeyError:
                            try:
                                res_structure = chain[res.resid]
                            except KeyError:
                                continue
                        if res_structure.resname != res.resname.upper():
                            continue
                # If we have a modified residue
                else:
                    for _res in chain:
                        if _res.get_id()[1] == res.auth_resid:
                            res_structure = _res
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
        # Set a residue as reference
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
    def _get_nearest_equivalent(self, other, reslist, site):
        equivalents = []
        for res in reslist:
            if res.structure is None:
                continue
            if res.is_equivalent(self):
                equivalents.append(res)
        result = None
        min_dist = 999
        for eq in equivalents:
            #Check if the same residue is already in the site
            if site.contains_equivalent(eq):
                continue
            dist = eq.get_distance(other, minimum=True)
            if dist<min_dist:
                result = eq
                min_dist = dist
        return result
    
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
            rms: RMSD after fitting, not including outliers
            rms_all: RMSD over all atoms, including outliers
        """
        min_rms = 999
        result_rot, result_tran, result_rms, result_rms_all = None, None, None, None
        p_all = np.array(p_coords, copy=True)
        q_all = np.array(q_coords, copy=True)
        # Initialize Biopython SVDSuperimposer
        sup = SVDSuperimposer()
        for i in range(cycles):
            if p_coords.size == 0:
                break
            sup.set(p_coords, q_coords)
            sup.run()
            rms = sup.get_rms()
            rot, tran = sup.get_rotran()
            if rms < min_rms:
                result_rot, result_tran, result_rms = rot, tran, rms
            # Transform coordinates
            q_trans = np.dot(q_coords, rot) + tran
            # Find outliers
            diff = np.linalg.norm(p_coords - q_trans, axis=1)
            to_keep = np.where(diff < cutoff)
            # Reject outliers
            p_coords = p_coords[to_keep]
            q_coords = q_coords[to_keep]
        # Also compute RMSD over all atoms
        q_all = np.dot(q_all, result_rot) + result_tran
        result_rms_all = PdbSite._rmsd(p_all, q_all)
        return result_rot, result_tran, np.round(result_rms, 3), np.round(result_rms_all, 3)

    @staticmethod
    def _rmsd(p_coords, q_coords):
        """Calculates rmsd on two coordinate sets (NumPy arrays) WITHOUT
        transformation and minimization"""
        diff = np.square(np.linalg.norm(p_coords - q_coords, axis=1))
        return np.sqrt(np.sum(diff) / diff.size)

    @staticmethod
    def _transform(coords, rot, tran):
        """Rotates and translates a set of coordinates (NxD NumPy array)"""
        return np.dot(coords, rot) + tran


class Box:

    def __init__(self, residue_list, headroom=1):
        """Define a box using the provided residues coordinates adding
        some safe distance around (headroom)"""
        self.residue_list = residue_list
        self.headroom = float(headroom)
        self.points = self._get_boundaries(self.headroom)

    def __contains__(self, residue):
        """Checks if the provided het (or residue) is in the box"""
        for atom in residue.structure.get_atoms():
            x,y,z = atom.get_coord()
            if self.points['min_x'] <= x <= self.points['max_x'] and \
                    self.points['min_y'] <= y <= self.points['max_y'] and \
                    self.points['min_z'] <= z <= self.points['max_z']:
                return True
        return False

    def mean_distance_from_residues(self, het):
        """Calculates the mean distance of an arbitrary residue (protein or HET)
        from the residues that define the box"""
        dist_sum = 0
        nofres = len(self.residue_list)
        for residue in self.residue_list:
            dist_sum += Box.min_distance(residue, het.structure)
        return round(dist_sum / nofres, 3)

    def similarity_with_cognate(self, het):
        """Checks the similarity score of the given compound with the cognate
        ligand of the given pdb, using the PARITY-derived data"""
        try:
            pdb_id = het.pdb_id
            hetcode = het.resname.upper()
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

