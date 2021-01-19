from .residue_definitions import AA_3TO1
from .UniResidue import UniResidue


class UniSite:
    """M-CSA UniProt catalytic site. Contains lists of UniResidues and mapped PDB
    catalytic sites (UniSite objects)."""

    def __init__(self):
        self.residues = []
        self.residues_dict = {}
        self.mapped_pdbsites = []
        self.reference_site = None

    def __str__(self):
        """Print as pseudo-sequence in one-letter code"""
        return self.sequence

    def __iter__(self):
        """Iterate over residues"""
        yield from self.residues

    def __contains__(self, residue):
        """Check if residue of specific ID is there"""
        return residue.full_id in self.residues_dict

    def __getitem__(self, full_id):
        """Return the child with given ID."""
        return self.residues_dict[full_id]

    # Alternative constructors

    @classmethod
    def build_reference(cls, reslist):
        """Builds reference active site from a list of UniProt catalytic residues."""
        ref = UniSite.from_list(reslist)
        ref.reference_site = ref
        return ref

    @classmethod
    def build(cls, reslist, reference_site=None):
        """Build UniProt site from a list of residues, and map reference residues"""
        site = cls.from_list(reslist)
        if site.is_reference:
            site.reference_site = site
        else:
            site.reference_site = reference_site
            site._map_reference_residues()
        return site

    @classmethod
    def from_list(cls, reslist):
        """Construct UniSite object directly from residue list"""
        site = cls()
        reslist = UniSite._cleanup_list(reslist)
        for res in reslist:
            site.add(res)
        return site

    # Properties

    @property
    def mcsa_id(self):
        """Get M-CSA ID of catalytic residues."""
        for res in self.residues:
            if res.mcsa_id:
                return res.mcsa_id
        return

    @property
    def uniprot_id(self):
        """UniProt site ID is the UniProt id of the sequence"""
        return self.residues[0].uniprot_id

    @property
    def id(self):
        """UniProt site ID is the UniProt id of the sequence"""
        return self.uniprot_id

    @property
    def sequence(self):
        """Show as pseudo-sequence in one-letter code"""
        return ''.join([AA_3TO1[res.resname] for res in self.residues])

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
                if 'main' in res.funclocs:
                    result = True
                    continue
            if not res.is_conserved and not res.is_conservative_mutation:
                return False
            if res.is_conservative_mutation:
                result = True
        return result

    def get_residues(self):
        """Iterate over residues"""
        yield from self.residues

    def get_gaps(self):
        """Returns an index of the gap positions (non-aligned residues)"""
        gaps = []
        for i, res in enumerate(self.residues):
            if res.is_gap:
                gaps.append(i)
        return gaps

    def add(self, residue):
        """Add UniProt residue to list"""
        if type(residue) == UniResidue:
            self.residues.append(residue)
            self.residues_dict[residue.full_id] = residue
        else:
            print('Attempted to add non-UniRes object in UniSite')
            return
        return True

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
                gap = UniResidue()
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
                if p.is_equivalent(q, by_chiral_id=False):
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
