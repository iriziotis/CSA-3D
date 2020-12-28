# !/usr/bin/env python3

from .residue_definitions import AA_3TO1
from .UniResidue import UniResidue


class UniSite:
    """M-CSA Uniprot catalytic site. A list of UniResidue objects"""

    def __init__(self):
        self.residues = []
        self.residues_dict = {}
        self.mapped_pdbsites = []
        self.reference_site = None

    def __str__(self):
        return ''.join([AA_3TO1[res.resname] for res in self.residues])

    def __iter__(self):
        """Iterate over residues"""
        yield from self.residues

    def __contains__(self, _id):
        """Check if residue of specific is there"""
        return id in self.residues_dict

    def __getitem__(self, _id):
        """Return the child with given id."""
        return self.residues_dict[id]

    def get_residues(self):
        yield from self.residues

    def add(self, residue):
        """Add UniProt residue to list"""
        if type(residue) == UniResidue:
            self.residues.append(residue)
            self.residues_dict[residue.id] = residue
        else:
            print('Attempted to add non-UniRes object in UniSite')
            return
        return True

    @classmethod
    def build(cls, reslist, reference_site=None):
        """Build UniProt site from a list of residues, and map reference residues"""
        site = cls.from_list(reslist)
        if site.is_reference:
            site.reference_site = None
        else:
            site.reference_site = reference_site
            site._map_reference_residues()
        return site

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

    @property
    def mcsa_id(self):
        """Get M-CSA ID of catalytic residues."""
        for res in self.residues:
            if res.mcsa_id:
                return res.mcsa_id
        return

    @property
    def id(self):
        """UniProt site ID is the UniProt id of the sequence"""
        return self.residues[0].uniprot_id

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

    @classmethod
    def from_list(cls, res_list):
        """Construct UniSite object directly from residue list"""
        site = cls()
        for res in res_list:
            site.add(res)
        return site
