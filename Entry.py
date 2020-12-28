# !/usr/bin/env python3

from .config import UNI2PDB
from .PdbSite import PdbSite
from .UniSite import UniSite


class Entry:
    """M-CSA entry. Contains a list of UniSites and a list of PdbSites"""

    def __init__(self, mcsa_id=None):
        self.mcsa_id = mcsa_id
        self.unisites = []
        self.unisites_dict = {}
        self.pdbsites = []
        self.pdbsites_dict = {}
        self.pdb_ids = set()

    def add(self, site):
        """Adds UniProt catalytic Site in unisites list"""
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
        """Return UniProt active site with id"""
        if _id in self.unisites_dict:
            return self.unisites_dict[_id]
        return

    def get_pdbsite(self, _id):
        """Return PDB active site with specific id"""
        if _id in self.unisites_dict:
            return self.unisites_dict[_id]
        return

    def get_unisites(self):
        """Return generator of all UniProt sites to iterate through"""
        yield from self.unisites

    def get_pdbsites(self, pdb_id=None):
        """Return PDB active sites that have the same PDB ID, or return a generator
        of all pdbsites to iterate through"""
        result = []
        if not pdb_id:
            yield from self.pdbsites
        else:
            result = [val for key, val in self.pdbsites_dict.items() if pdb_id in key]
        yield from result

    def update_mapped_sites(self):
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
