import os
import numpy as np
import pandas as pd
from copy import deepcopy
from collections import defaultdict
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import linkage, cut_tree
from .config import UNI2PDB
from .PdbSite import PdbSite
from .UniSite import UniSite


class Entry:
    """M-CSA entry. Contains a list of UniSite objects and a list
    of PdbSite objects"""

    def __init__(self, mcsa_id=None):
        self.mcsa_id = mcsa_id
        self.info = None
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


    def rmsd_matrix(self, sitelist):
        """Contstructs RMSD matrix for the sites provided in sitelist"""
        rmsds = []
        ids = []
        seen = set()
        for p in sitelist:
            ids.append(p.id)
            for q in sitelist:
                if p.id == q.id or (q.id, p.id) in seen:
                    continue
                seen.add((p.id, q.id))
                try:
                    _, _, _, rms_all = p.fit(q, weighted=True)
                except Exception:
                    _, _, _, rms_all = p.fit(q, weighted=False, cycles=10, cutoff=5)
                rmsds.append(rms_all)
        rmsds = np.array(rmsds, dtype='float32')
        matrix = pd.DataFrame(squareform(rmsds), columns=ids, index=ids)
        return matrix

    def clustering(self, matrix, outdir=None, height=None):
        """Performs hierarchical clustering on the provided RMSD matrix."""
        ids = list(matrix.index)
        # Compute linkage matrix
        Z = linkage(matrix, method='average')
        # Cut tree at automatic height
        if not height:
            height = 0.7 * max(Z[:, 2])
        # Get clusters index
        clusters = cut_tree(Z, height=height).flatten()
        # Group sites to clusters
        cluster_dict = defaultdict(list)
        for cluster in np.unique(clusters):
            for i, id in zip(clusters, ids):
                if i == cluster:
                    cluster_dict[cluster].append(self.get_pdbsite(id))
        return Z, cluster_dict

    def create_template(self, ca=False, subset=None):
        """Creates template from conserved sites"""

        # Get reference as functional site (maybe make a separate method for retrieving reference from entry)
        for site in self.get_pdbsites():
            if site.is_reference:
                reference = site.get_functional_site(ca=ca)

        # Calculate average coordinates
        avg = reference.copy()
        for site in self.get_pdbsites(sane_only=True):
            if subset and site.id not in subset:
                continue
            if not (site.is_conserved or site.is_conservative_mutation) or site.is_reference:
                continue
            funcsite = site.get_functional_site(ca=ca)
            rot, tran, rms, rms_all = avg.fit(funcsite, transform=True, ca=ca)
            site_coords = np.array([atom.get_coord() for res in funcsite for atom in res.structure])
            avg_coords = np.array([atom.get_coord() for res in avg for atom in res.structure])
            avg_coords = np.mean([avg_coords, site_coords], axis=0)
            for atom, coord in zip([atom for res in avg for atom in res.structure], avg_coords):
                atom.set_coord(coord)    

        # Find representative site
        min_rms = 999
        template = None
        for site in self.get_pdbsites(sane_only=True):
            if not (site.is_conserved or site.is_conservative_mutation) or site.is_reference:
                continue
            rot, tran, rms, rms_all = avg.fit(site, transform=True, ca=ca)
            site.write_pdb(func_atoms_only=True, write_hets=False)
            if rms_all < min_rms:
                min_rms = rms_all
                template = site

        # Fit template to reference
        reference.fit(template, transform=True)
        reference.fit(avg, transform=True)

        return template, avg


