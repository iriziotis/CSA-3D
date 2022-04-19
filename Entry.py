import os
import numpy as np
import pandas as pd
from copy import deepcopy
from collections import defaultdict
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import linkage, cut_tree
from dendrogram_cut import DendrogramCut
from .config import UNI2PDB
from .residue_definitions import RESIDUE_DEFINITIONS, TEMPLATE_FUZZY_RESIDUES, TEMPLATE_FUZZY_ATOMS, AA_3TO1
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

    def clustering(self, matrix, height=None, get_linkage_matrix=False):
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
                    cluster_dict[cluster].append(self.get_pdbsite(id).id)
        if get_linkage_matrix:
            return cluster_dict, Z
        return cluster_dict

    def clustering_bayesian(self, matrix, k_max=4, l=1.0, plot_outfile=None):
        """Performs hierarchical clustering on the provided RMSD matrix,
        and cuts the tree automatically using Bayesian statistics"""
        cluster_dict = defaultdict(list)
        ids = list(matrix.index)
        matrix = np.array(matrix)
        model = DendrogramCut(k_max=k_max, method='average').fit(matrix)
        k = model.pac_bayesian_cut(lambda_=l)
        clusters = model.get_cluster_label(k=k)
        if plot_outfile:
            fig = model.heatmap_with_dendrogram_plot(k=k)
            fig.write_image(plot_outfile)
        for i,cluster_no in enumerate(clusters):
            cluster_dict[cluster_no].append(ids[i])
        return cluster_dict

    def create_template(self, ca=False, outdir=None, outfile=None, subset=None, cluster_no=None):
        """Creates template from conserved sites"""
        # Get reference as functional site (maybe make a separate method for retrieving reference from entry)
        for site in self.get_pdbsites():
            if site.is_reference:
                reference = site.get_functional_site(ca=ca)
                break
        # Calculate average coordinates
        avg = reference.copy()
        for site in self.get_pdbsites(sane_only=True):
            if subset and site.id not in subset:
                continue
            if not (site.is_conserved or site.is_conservative_mutation) or site.is_reference:
                continue
            funcsite = site.get_functional_site(ca=ca)
            avg.fit(funcsite, transform=True, ca=ca)
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
            rot, tran, rms, rms_all = avg.fit(site, transform=False, ca=ca)
            if rms_all < min_rms:
                min_rms = rms_all
                template = site
        # Fit template to reference
        reference.fit(template, transform=True)
        # Write template
        self.write_template(template, subset, cluster_no, outdir, outfile)
        return template, avg

    def write_template(self, template, subset=None, cluster_no=None, outdir=None, outfile=None):
        """
        Writes template coordinates and constraints in TESS/Jess format
        Args:
            outdir: Output directory
            outfile: Default 'mcsa_id.cluster_no.template.pdb' 
        """
        if not outdir:
            outdir = '.'
        if not outfile:
            cluster = 'all' if cluster_no is None else f'cluster_{cluster_no}'
            outfile = f'{outdir}/csa3d_{str(template.mcsa_id).zfill(4)}.{cluster}.template.pdb'
        with open(outfile, 'w') as o:
            remarks = (f'REMARK TEMPLATE\n'
                       f'REMARK MCSA_ID {template.mcsa_id}\n'
                       f'REMARK PDB_ID {template.pdb_id}\n'
                       f'REMARK UNIPROT_ID {template.uniprot_id}\n'
                       f'REMARK EC {template.ec}\n'
                       f'REMARK ENZYME {template.enzyme}\n'
                       f'REMARK EXPERIMENTAL_METHOD {template.experimental_method}\n'
                       f'REMARK RESOLUTION {template.resolution}\n'
                       f'REMARK ORGANISM_NAME {template.organism_name}\n'
                       f'REMARK ORGANISM_ID {template.organism_id}\n')
            print(remarks, file=o)
            alt_residues = self.get_alt_residues(template)
            matchcodes = self.get_matchnumbers(template, alt_residues)
            dist_cutoffs = self.get_dist_cutoffs(template, subset=subset)
            for i, res in enumerate(template):
                if res.structure is not None:
                    if res.has_main_chain_function or not res.is_standard:
                        resname = 'ANY'
                    else:
                        resname = res.structure.get_resname()
                    for atom in res.structure.get_atoms():
                        funcstring = '{}.{}'.format(resname, atom.get_id().upper())
                        if funcstring not in RESIDUE_DEFINITIONS:
                            continue
                        pdb_line = '{:6}{:5d} {:<4}{}{:>3}{:>2}{:>4}{:>12.3f}{:>8.3f}{:>8.3f} {:<4.4}{:<5.2f}'.format(
                            'ATOM', matchcodes[int(i)][atom.name], 
                            atom.name if len(atom.name) == 4 else ' {}'.format(atom.name), 'Z',
                            resname, res.structure.get_parent().get_id(), res.structure.get_id()[1],
                            atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2],
                            alt_residues[i], dist_cutoffs[i])
                        print(pdb_line, file=o)
            print('END', file=o)

    def get_alt_residues(self, template):
        """Returns a string of alternative residues to be matched in the template"""
        observed = defaultdict(set)
        expected = defaultdict(set)
        intersection = defaultdict(set)
        # Get alternative residues in each position from sequence alignment
        for site in template.parent_entry.unisites:
            if site.is_conserved or site.is_conservative_mutation:
                for i, res in enumerate(site):
                    observed[i].add(AA_3TO1[res.resname])
        # Compare expected with observed substitutions
        ref = template.reference_site
        for i, res in enumerate(ref):
            expected[i] = set(list(TEMPLATE_FUZZY_RESIDUES.get(res.resname.upper(), '')))
            expected[i].add(AA_3TO1[res.resname])
        for i in expected:
            intersection[i] = ''.join(list(expected[i].intersection(observed[i])))
        return intersection

    def get_matchnumbers(self, template, fuzzy_residues):
        """Returns the Jess match option number for each template atom"""
        matchnumbers = defaultdict(dict)
        for i, res in enumerate(template):
            resname = res.resname.upper()
            if 'ptm' in res.funclocs:
                resname = 'PTM'
            for funcloc in res.funclocs:
                if 'main' in funcloc:
                    resname = 'ANY'
            for atom in res.structure:
                if atom.name in TEMPLATE_FUZZY_ATOMS[resname]:
                    matchnumbers[i][atom.name] = TEMPLATE_FUZZY_ATOMS[resname][atom.name][0]
        return matchnumbers

    def get_dist_cutoffs(self, template, subset=None, 
            comparisons_file='/Users/riziotis/ebi/phd/datasets/csa3d/variation/data/per_entry/csa3d_0133.csv'):
        """Returns the distance threshold weight for each residue"""
        comparisons = pd.read_csv(open(comparisons_file, 'r'))
        if subset:
            comparisons = comparisons.query('p_id in @subset and q_id in @subset')
        per_res = comparisons.loc[~pd.isnull(comparisons['per_res_rms']), 'per_res_rms'].astype('str').str.split(';', expand=True).astype('float32')
        cutoffs = dict(per_res.describe().loc[['mean', 'std']].sum())
        return cutoffs

            
        








