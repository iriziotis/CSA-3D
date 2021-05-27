import warnings
import numpy as np
import Bio.PDB.Selection
from copy import copy
from collections import defaultdict
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.MMCIFParser import MMCIFParser, FastMMCIFParser
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
from Bio.PDB.NeighborSearch import NeighborSearch
from rmsd import reorder_hungarian
from .residue_definitions import AA_3TO1, RESIDUE_DEFINITIONS, EQUIVALENT_ATOMS
from .config import PDB2EC, PDB2UNI
from .PdbResidue import PdbResidue, Het
from .Superimposer import Superimposer

warnings.simplefilter('ignore', PDBConstructionWarning)

class PdbSite:
    """M-CSA PDB catalytic site. Contains lists of PdbResidues and mapped UniProt
    catalytic sites (UniSite objects), a 3D structure (Biopython Structure) built
    from individual PdbResidue structures (Biopython Residue), a parent structure
    (Biopython Structure), all and close to the site hetero components (possible
    ligands according to their chemical similarity to the cognate ligand and their
    centrality in the active site) as Het objects (containing a Biopython Residue
    structure), as well as a dictionary of annotations extracted from the parent
    mmCIF assembly structure file and SIFTS"""

    def __init__(self):
        self.parent_entry = None
        self.residues = []
        self.residues_dict = {}
        self.mapped_unisites = []
        self.reference_site = None
        self.parent_structure = None
        self.structure = None
        self.ligands = []
        self.mmcif_dict = dict()
        self.is_sane = None

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
                if res.full_id not in self.residues_dict:
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
    def from_list(cls, reslist, cif_path, parent_entry, annotate=True):
        """Construct PdbSite object directly from residue list"""
        mmcif_dict = dict()
        # First reduce redundant residues with multiple function locations
        reslist = PdbSite._cleanup_list(reslist)
        site = cls()
        site.parent_entry = parent_entry
        try:
            if annotate: 
                parser = MMCIFParser(QUIET=True) 
                structure = parser.get_structure('', cif_path)
                mmcif_dict = parser._mmcif_dict
            else:
                parser = FastMMCIFParser(QUIET=True)
                structure = parser.get_structure('', cif_path)
        except (TypeError, PDBConstructionException):
            warnings.warn('Could not build site from residue list. Check entry', RuntimeWarning)
            return
        for res in reslist:
            if structure:
                res.add_structure(structure)
            site.add(res)
        if annotate:
            site.parent_structure = structure
            site.mmcif_dict = mmcif_dict
            site.find_ligands()
        return site

    @classmethod
    def build_reference(cls, reslist, parent_entry, cif_path, annotate=True):
        """Builds reference active site from a list of PDB catalytic residues.
        Assumes that the list only contains one active site, so use it only
        if it is a list of manually annotated catalytic residues"""
        ref = PdbSite.from_list(reslist, cif_path, parent_entry, annotate)
        ref.reference_site = ref
        ref.is_sane = True
        return ref

    @classmethod
    def build(cls, seed, reslist, reference_site, parent_entry):
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
        site.parent_entry = parent_entry
        site._map_reference_residues()
        return site

    @classmethod
    def build_all(cls, reslist, reference_site, parent_entry, cif_path, annotate=True, redundancy_cutoff=None):
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
        except (TypeError, PDBConstructionException):
            warnings.warn('Could not parse structure {}'.format(cif_path, RuntimeWarning))
            return sites
        # First reduce redundant residues with multiple function locations
        reslist = PdbSite._cleanup_list(reslist)
        # We want all equivalent residues from identical assembly chains
        reslist = PdbSite._get_assembly_residues(reslist, structure)
        # Get seeds to build active sites
        seeds = PdbSite._get_seeds(reslist)
        # Build a site from each seed
        for seed in seeds:
            site = cls.build(seed, reslist, reference_site, parent_entry)
            sites.append(site)
            if site.has_missing_functional_atoms or len(site) != len(reference_site):
                sites.remove(site)
                continue
            # Reduce redundancy
            try:
                _, _, _, rms_all = site.fit(sites[-2])
                if site.has_identical_residues(sites[-2]) and rms_all == 0:
                    sites.remove(site)
                    continue
                if redundancy_cutoff:
                    if rms_all < redundancy_cutoff:
                        sites.remove(site)
                        continue
            except IndexError:
                pass
            # Add ligands and annotations
            if annotate and structure:
                site.parent_structure = structure
                site.mmcif_dict = mmcif_dict
                site.find_ligands()
        # Final clean up of unclustered sites
        PdbSite._mark_unclustered(sites)
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
                    try:
                        return PDB2EC[(self.pdb_id, res.alt_chain[0])]
                    except KeyError:
                        continue
        return

    @property
    def sequence(self):
        """Show as pseudo-sequence in one-letter code"""
        return ''.join([AA_3TO1[res.resname] if (res.is_standard or res.is_gap) else 'X' for res in self])

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
    def acts_on_polymer(self):
        """Check if it belongs to a family of enzymes whose substrate is a polymer
        (protein or nucleic)"""
        return self.parent_entry.info['reaction']['is_polymeric']

    @property
    def is_reference(self):
        """Check if site is reference site"""
        if self.size > 0:
            return self.residues[0].is_reference
        return False

    @property
    def is_conserved(self):
        """Check if all residues are conserved by comparing to the reference"""
        if self.is_reference:
            return True
        return str(self) == str(self.reference_site)

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
        if self.is_reference:
            return all([i.structure for i in ref])
        try:
            gaps = set(self.get_gaps())
            self_atoms, _ = self._get_func_atoms(omit=gaps)
            ref_atoms, _ = self.reference_site._get_func_atoms(omit=gaps)
            return len(self_atoms) != len(ref_atoms)
        except (TypeError, ValueError):
            return True

    # Methods

    def copy(self, include_structure=True):
        """Returns a copy of the site. If include_structure is False,
        then the structure is not copied"""
        site = copy(self)
        if include_structure:
            site.structure = self.structure.copy()
        return site

    def add(self, residue):
        """Add PdbResidue object to site (in the residues list and dict)"""
        residue = residue.copy(include_structure=True)
        if type(residue) == PdbResidue:
            self.residues.append(residue)
            self.residues_dict[residue.full_id] = residue
            residue.parent_site = self
        if type(residue) == Het:
            self.ligands.append(residue)
            residue.parent_site = self
            if residue.is_polymer:
                if residue.chain in self.structure[0]:
                    for r in residue.structure:
                        self.structure[0][residue.chain].add(r)
                    return True
                self.structure[0].add(residue.structure)
                return True
        if residue.structure:
            # Initialize structure if empty
            if self.structure is None:
                self.structure = Structure(self.id)
                self.structure.add(Model(0))
            chain_id = residue.structure.get_parent().get_id()
            if chain_id not in self.structure[0]:
                self.structure[0].add(Chain(chain_id))
            # Add residue structure to site structure
            if residue.structure.get_id() not in self.structure[0][chain_id]:
                self.structure[0][chain_id].add(residue.structure)
        return True

    def get_distances(self, kind='com'):
        """Calculates all intra-site residue distances and returns a
        numpy array"""
        dists = []
        seen = set()
        for p in self.residues:
            for q in self.residues:
                if p == q or (q.full_id, p.full_id) in seen:
                    continue
                if p.is_gap or q.is_gap:
                    dists.append(np.nan)
                else:   
                    dists.append(p.get_distance(q, kind))
                seen.add((p.full_id, q.full_id))
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
        for i, p in enumerate(self):
            for j, q in enumerate(self):
                if p==q or q.is_gap or q.is_gap:
                    continue
                if p.is_equivalent(q, by_chiral_id=False, by_chain=False):
                    if (j, i) not in identicals:
                        identicals.add((i, j))
        return identicals
    

    def find_ligands(self, radius=3):
        """
        Searches the parent structure for hetero components close to the
        catalytic residues, by searching around the atoms of catalytic residues
        and the dummy atoms between distant residues. Populates the ligands list 
        with Het objects.

        Args:
            radius: the search space (in Å) around the atoms of the catalytic residues
        """
        if type(self.parent_structure) != Structure:
            return
        # Get centers of search
        centers = self._get_ligand_search_centers(radius)        
        # Initialize KD tree
        query_atoms = Bio.PDB.Selection.unfold_entities(self.parent_structure, 'A')
        ns = NeighborSearch(query_atoms)
        # Search for ligands around each center
        polymers = defaultdict(list)
        site_chains = set([res.chain for res in self])
        seen = set()
        for center in centers:
            hits = ns.search(center, radius, level='R')
            for res in hits:
                if res.get_full_id() in seen:
                    continue
                seen.add(res.get_full_id())
                restype = res.get_id()[0][0]
                chain = res.get_parent().get_id()
                # Ignore waters
                if restype == 'W':
                    continue
                # HET components
                if restype == 'H':
                    self.add(Het(self.mcsa_id, self.pdb_id, res.get_resname(), 
                                 res.get_id()[1], chain, structure=res, parent_site=self))
                # Protein/nucleic polymer components
                if restype == ' ' and chain not in site_chains:
                    polymers[chain].append(res)
        # Build polymers
        if self.acts_on_polymer:
            for chain, reslist in polymers.items():
                self.add(Het.polymer(reslist, self.mcsa_id, self.pdb_id, chain, self))
        return 

    def write_pdb(self, outdir=None, outfile=None, write_hets=False, func_atoms_only=False, include_dummy_atoms=False):
        """
        Writes site coordinates in PDB format
        Args:
            write_hets: Include coordinates of ligands.
            outdir: Directory to save the .pdb file
            outfile: If unspecified, name is formatted to include info on M-CSA ID, 
                     chain of each catalytic residue, annotation if the site is a
                     reference site and an annotation about the conservation, relatively
                     to the reference (c: conserved, m: mutated, cm: has only conservative
                     mutation)
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
            if self.is_sane:
                sanity = 'sane'
            else:
                sanity = 'insane'
            outfile = '{}/mcsa_{}.{}.{}.{}.{}.{}.pdb'.format(
                      outdir.rstrip('/'), str(self.mcsa_id).zfill(4), self.id,
                      'reference' if self.is_reference else 'cat_site', 
                      conservation, atms, sanity)
        with open(outfile, 'w') as o:
            if bool(self.mmcif_dict):
                ligands = ','.join('{0.resname};{0.resid};{0.chain};{0.similarity};{0.centrality};{0.type}'.format(h) for h in self.ligands)
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
                           'REMARK NEARBY_LIGANDS {1}'.format(self, ligands))
                print(remarks, file=o)
            residues = self.residues.copy()
            if write_hets:
                residues += self.ligands
            for res in residues:
                if not include_dummy_atoms and res.is_gap:
                    continue
                structure = res.structure
                if res.dummy_structure:
                    structure = res.dummy_structure
                if structure is not None:
                    for atom in structure.get_atoms():
                        resname = res.resname.upper()
                        if res.has_main_chain_function or not res.is_standard:
                            resname = 'ANY'
                        funcstring = '{}.{}'.format(resname, atom.get_id().upper())
                        if func_atoms_only and type(res) == PdbResidue and funcstring not in RESIDUE_DEFINITIONS:
                            continue
                        pdb_line = '{:6}{:5d} {:<4}{}{:>3}{:>2}{:>4}{:>12.3f}' \
                                   '{:>8.3f}{:>8.3f} {:6}'.format(
                            'HETATM' if (atom.get_parent().get_id()[0] != ' ' or type(res) == Het) else 'ATOM',
                            int(atom.get_serial_number()) if atom.get_serial_number() else 0,
                            atom.name if len(atom.name) == 4 else ' {}'.format(atom.name),
                            'Z' if funcstring in RESIDUE_DEFINITIONS else atom.get_altloc(),
                            atom.get_parent().get_resname(),
                            atom.get_parent().get_parent().get_id(),
                            atom.get_parent().get_id()[1],
                            atom.get_coord()[0],
                            atom.get_coord()[1],
                            atom.get_coord()[2],
                            atom.get_occupancy() if atom.get_occupancy() else '')
                        print(pdb_line, file=o)
            print('END', file=o)

    def fit(self, other, weighted=False, cycles=1, cutoff=999, scaling_factor=None, transform=False, 
            mutate=True, reorder=True, allow_symmetrics=True, exclude=None, get_array=False):
        """Iteratively fits two catalytic sites (self: fixed site, other: mobile site)
        using the Kabsch algorithm from the rmsd module (https://github.com/charnley/rmsd).
        Can also find the optimal atom alignment in each residue, considering
        symmetrical atoms and functionally similar residues, using the
        Hungarian algorithm.

        Args:
            other: mobile active site to fit
            weighted: to perform weighted superposition in the last iteration
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
        # If we want to exclude residues from fitting
        if exclude is not None:
            if type(exclude) not in (list, tuple, set):
                exclude = [exclude]
            for i in exclude:
                gaps.add(i)
        # Get atom identifier strings and coords as numpy arrays
        p_atoms, p_coords = self._get_func_atoms(allow_symmetrics, omit=gaps)
        q_atoms, q_coords = other._get_func_atoms(allow_symmetrics, omit=gaps)
        if p_atoms is None or q_atoms is None:
            return None, None, None, None
        if len(p_atoms) != len(q_atoms):
            raise Exception('Atom number mismatch in sites {} and {}'.format(self.id, other.id))
        # Initial crude superposition
        rot, tran, rms, _ = PdbSite._super(p_coords, q_coords, cycles=1)
        q_trans = PdbSite._transform(q_coords, rot, tran)
        # In case of non-conservative mutations, make pseudo-mutations to facilitate superposition
        if mutate:
            for i, (p_atom, q_atom) in enumerate(zip(p_atoms, q_atoms)):
                if p_atom != q_atom:
                    #q_atoms[i] = p_atom 
                    q_atoms[i] = '{}.MUT'.format(q_atoms[i].split('.')[0])
                    p_atoms[i] = '{}.MUT'.format(p_atoms[i].split('.')[0])
        # Reorder atoms using the Hungarian algorithm from rmsd package
        if reorder:
            q_review = reorder_hungarian(p_atoms, q_atoms, p_coords, q_trans)
            q_coords = q_coords[q_review]
        # Iterative superposition. Get rotation matrix, translation vector and RMSD
        rot, tran, rms, rms_all = PdbSite._super(p_coords, q_coords, cycles, cutoff, weighted, scaling_factor)
        if transform:
            other.structure.transform(rot, tran)
        if get_array:
            q_trans = np.dot(q_coords, rot) + tran
            return rot, tran, rms, rms_all, p_coords, q_trans
        return rot, tran, rms, rms_all

    def per_residue_rms(self, other, rot=None, tran=None, transform=False):
        """Calculates the RMSD of each residue in two superimposed sites.
        If superposition rotation matrix and translation vector are not given,
        RMSD is calculated without transformation. Otherwise, fitting is performed
        automatically, using weighted superposition to compensate for bias caused
        by slightly outlying residues."""
        rmsds = []
        if rot is None or tran is None:
            rot, tran, _, _ = self.fit(other, weighted=True, transform=False)
        for i, (p, q) in enumerate(zip(self, other)):
            if p.is_gap or q.is_gap:
                rmsds.append(np.nan)
                continue
            # Get functional atoms
            p_atoms, p_coords = p.get_func_atoms()  
            q_atoms, q_coords = q.get_func_atoms()  
            # Mutate if there are mismatches
            for i, (p_atom, q_atom) in enumerate(zip(p_atoms, q_atoms)):
                if p_atom != q_atom:
                    p_atoms[i] = 'MUT'
                    q_atoms[i] = 'MUT'
            # Transform functional atoms
            if transform:
                q_coords = PdbSite._transform(q_coords, rot, tran)
            # Reorder
            q_review = reorder_hungarian(p_atoms, q_atoms, p_coords, q_coords)
            q_coords = q_coords[q_review]
            # Calculate RMSD
            rms = PdbSite._rmsd(p_coords, q_coords)
            rmsds.append(np.round(rms, 3))
        return np.array(rmsds)

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
                gap = PdbResidue(mcsa_id=self.mcsa_id, 
                                 pdb_id=self.pdb_id, 
                                 chiral_id=reference_residue.chiral_id,
                                 dummy_structure=True)
                gap.reference_residue = reference_residue
                self.add(gap)
        self._reorder()
        return

    def _get_func_atoms(self, allow_symmetrics=True, omit=None):
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
                return np.array(atoms), np.array(coords)
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
        try:
            atoms = np.array(atoms, dtype=object)
            coords = np.stack(coords, axis=0)
        except ValueError:
            return None, None
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
    
    def _get_ligand_search_centers(self, radius=4):
        """Gets atom coordinates from catalytic residues, and interpolates the empty 
        space between distant residues, by calculating the center of geometry of the
        two residues. Radius is used to identify distant residues in between which an
        extra center will be added."""
        centers = []
        seen = set()
        for p in self:
            if p.is_gap or p.structure is None:
                continue
            for atom in p.structure.get_unpacked_list():
                centers.append(atom.get_coord()) 
            p_centroid = p.structure.center_of_mass(geometric=True)
            for q in self:
                if p is q or (q.id, p.id) in seen or q.is_gap or q.structure is None:
                    continue
                seen.add((p.id, q.id))
                q_centroid  = q.structure.center_of_mass(geometric=True)
                dist = p.get_distance(q, kind='min')
                if 2*radius <= dist <= 4*radius:
                    dummy_coords = np.mean([p_centroid, q_centroid], axis=0)
                    centers.append(dummy_coords)
        return centers

    @staticmethod
    def _cleanup_list(reslist):
        """Finds duplicate residues of different funclocs and makes a single
        one with two funclocs. Returns a new list without redundant residues"""
        new_reslist = []
        seen = set()
        ignore = set()
        for p in reslist:
            for q in reslist:
                if p==q or (q.full_id, p.full_id) in seen:
                    continue
                if p.is_equivalent(q, by_chiral_id=False, by_chain=True):
                    if p.funclocs != q.funclocs:
                        new_res = p.copy(include_structure=True)
                        new_res.funclocs = [p.funclocs[0], q.funclocs[0]]
                        new_reslist.append(new_res)
                        ignore.add(p.full_id)
                        ignore.add(q.full_id)
                seen.add((p.full_id, q.full_id))
        for p in reslist:
            if p.full_id not in ignore and p not in new_reslist:
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
                new_res = res.copy(include_structure=False)
                new_res.chain = chain.get_id()
                new_res.structure = res_structure
                new_reslist.append(new_res)
        return new_reslist

    @staticmethod
    def _get_seeds(reslist):
        """Finds residues in a list that can be used as seeds when
        building multiple active sites"""
        seeds = []
        # Set a residue as reference
        ref = None
        for res in reslist:
            if res.auth_resid is None or res.structure is None:
                continue
            # Check if residue has any close neighbours -- If not, skip it
            skip = True
            for other in reslist:
                if res == other:
                    continue
                try:
                    if res.get_distance(other, kind='min') < 8:
                        skip = False
                        break
                except TypeError:
                    continue
            if skip:
                continue
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
        """Gets the closest equivalent of 'other' to 'self', if there
        are multiple equivalents in the residue list"""
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
            dist = eq.get_distance(other, kind='min')
            if dist<min_dist:
                result = eq
                min_dist = dist
        return result
    
    @staticmethod
    def _mark_unclustered(sitelist):
        """Returns a clean list of catalytic sites from the same PDB
        by rejecting sites that might have insanely outlying residues"""
        try:
            ref_dists = sitelist[0].reference_site.get_distances(kind='min')
            ref_dists = np.nan_to_num(ref_dists, nan=999)
            ref_dists = np.where(ref_dists < 8, 8, ref_dists)
        except IndexError:
            return False
        for p in sitelist:
            p.is_sane = True
            p_dists = np.nan_to_num(p.get_distances(kind='ca'), nan=0)
            if not np.all((p_dists < 3*ref_dists)):
                p.is_sane = False
                continue
            else:
                for q in sitelist:
                    if p.id == q.id or q.is_sane == False:
                        continue
                    q_dists = np.nan_to_num(q.get_distances(kind='ca'), nan=999)
                    q_dists = np.where(q_dists < 8, 8, q_dists)
                    if not np.all((p_dists < 1.3*q_dists)):
                        p.is_sane = False
        return True

    @staticmethod
    def _super(p_coords, q_coords, cycles=1, cutoff=999, weighted=False, scaling_factor=None):
        sup = Superimposer()
        sup.set(p_coords, q_coords, cycles, cutoff, scaling_factor)
        if weighted:
            sup.run_weighted()
        else:
            sup.run_unweighted()
        return sup.rot, sup.tran, np.round(sup.rms, 3), np.round(sup.rms_all, 3)

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

