from .config import UNI2PDB
from .PdbSite import PdbSite
from .UniSite import UniSite


class Entry:
    """M-CSA entry. Contains a list of UniSite objects and a list
    of PdbSite objects"""

    def __init__(self, mcsa_id=None):
        self.mcsa_id = mcsa_id
        self.unisites = []
        self.unisites_dict = {}
        self.pdbsites = []
        self.pdbsites_dict = {}
        self.pdb_ids = set()

    @property
    def has_multiple_pdb_refs(self):
        """Checks if entry has more than one PDB reference structures"""
        refs = set()
        for site in self.pdbsites:
            if site.is_reference:
                refs.add(site.pdb_id)
        return len(refs)>1

    @property
    def has_multiple_uniprot_refs(self):
        """Checks if entry has more than one PDB reference structures"""
        refs = set()
        for site in self.unisites:
            if site.is_reference:
                refs.add(site.uniprot_id)
        return len(refs)>1

    def add(self, site):
        """Adds catalytic site in the appropriate list depending on its
        type (PdbSite or UniSite)"""
        if type(site) == UniSite:
            self.unisites.append(site)
            self.unisites_dict[site.id] = site
        elif type(site) == PdbSite:
            self.pdbsites.append(site)
            self.pdbsites_dict[site.id] = site
            self.pdb_ids.add(site.pdb_id)
        else:
            print('Attempted to add non-Site object in entry')
            return
        self.update_mapped_sites()
        return True

    def get_unisite(self, _id):
        """Return UniProt active site with the given UniProt accession ID"""
        if _id in self.unisites_dict:
            return self.unisites_dict[_id]
        return

    def get_pdbsite(self, _id):
        """Return PDB active site with specific ID"""
        if _id in self.pdbsites_dict:
            return self.pdbsites_dict[_id]
        return

    def get_unisites(self):
        """Return a generator of all UniProt sites"""
        yield from self.unisites

    def get_pdbsites(self, pdb_id=None, sane_only=False):
        """Return a generator of PDB active sites that have the same PDB ID,
        or a generator of all PdbSites"""
        result = []
        if not pdb_id:
            for site in self.pdbsites:
                if sane_only and not site.is_sane:
                    continue
                yield site
        else:
            result = [val for key, val in self.pdbsites_dict.items() if pdb_id in key]
        for site in result:
            if sane_only and not site.is_sane:
                continue
            yield site

    def update_mapped_sites(self):
        """Updates the mapping of all UniProt catalytic sites to their
        respective PDB sites"""
        for unisite in self.unisites:
            if unisite.id in UNI2PDB:
                for sifts_pdb_id in UNI2PDB[unisite.id]:
                    if sifts_pdb_id in self.pdb_ids:
                        pdbsites = self.get_pdbsites(sifts_pdb_id)
                        for pdbsite in pdbsites:
                            if unisite not in pdbsite.mapped_unisites:
                                pdbsite.mapped_unisites.append(unisite)
                            if pdbsite not in unisite.mapped_pdbsites:
                                unisite.mapped_pdbsites.append(pdbsite)
        return True

    def all_vs_all(self, sane_only=False):
        """Returns a generator with all combinations of PDB sites in the entry.
        If sane_only, insane sites are excuded."""
        seen = set()
        for p in self.pdbsites:
            if sane_only:
                if not p.is_sane:
                    continue
            for q in self.pdbsites:
                if sane_only:
                    if not q.is_sane:
                        continue
                if p==q or (q.id, p.id) in seen:
                    continue
                seen.add((p.id, q.id))
                yield p, q
