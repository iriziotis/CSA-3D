import os
import json
import numpy as np
import pandas as pd
from copy import deepcopy
from collections import defaultdict
from itertools import combinations
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import linkage, cut_tree
from sklearn.cluster import KMeans
from dendrogram_cut import DendrogramCut
from .config import UNI2PDB, REF_RES_INFO
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

    def rmsd_matrix(self, sitelist, ca=False):
        """Contstructs RMSD matrix for the sites provided in sitelist"""
        rmsds = []
        ids = []
        seen = set()
        for p in sitelist:
            if not ca and p.has_missing_functional_atoms:
                continue
            ids.append(p.id)
            for q in sitelist:
                if not ca and q.has_missing_functional_atoms:
                    continue
                if p.id == q.id or (q.id, p.id) in seen:
                    continue
                seen.add((p.id, q.id))
                try:
                    _, _, _, rms_all = p.fit(q, weighted=True, ca=ca)
                except Exception:
                    _, _, _, rms_all = p.fit(q, weighted=False, ca=ca, cycles=10, cutoff=5)
                rmsds.append(rms_all)
        rmsds = np.array(rmsds, dtype='float32')
        matrix = pd.DataFrame(squareform(rmsds), columns=ids, index=ids)
        return matrix

    def clustering(self, matrix, height=None, get_linkage_matrix=False):
        """Performs hierarchical clustering on the provided RMSD matrix."""
        ids = list(matrix.index)
        if len(ids)<2:
            return {0: ids}
        # Compute linkage matrix
        Z = linkage(squareform(matrix), method='average')
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

    def clustering_bayesian(self, matrix, k_max=5, l=1.0, plot_outfile=None):
        """Performs hierarchical clustering on the provided RMSD matrix,
        and cuts the tree automatically using Bayesian statistics"""
        cluster_dict = defaultdict(list)
        ids = list(matrix.index)
        matrix = np.array(matrix)
        try:
            model = DendrogramCut(k_max=k_max, method='average').fit(matrix)
        except ValueError:
            return {0: ids}
        try:
            k = model.pac_bayesian_cut(lambda_=l)
        except RecursionError:
            return {0: ids}
        clusters = model.get_cluster_label(k=k)
        if plot_outfile:
            fig = model.heatmap_with_dendrogram_plot(k=k)
            fig.write_image(plot_outfile)
        for i,cluster_no in enumerate(clusters):
            cluster_dict[cluster_no].append(ids[i])
        return cluster_dict

    def _get_template_reference(self, ca=False):
        """Gets a reference site on which templates will be fitted into"""
        ref = None
        for site in self.get_pdbsites():
            if site.is_reference: 
                if ca or (not ca and not site.has_missing_functional_atoms):
                    ref = site
                    break
        if not ref:
            for site in self.get_pdbsites(sane_only=True):
                if (site.is_conserved or site.is_conservative_mutation): 
                    if ca or (not ca and not site.has_missing_functional_atoms):
                        ref = site
                        break
                    if not site.has_missing_functional_atoms:
                        ref = site
                        break
        return ref

    def create_template(self, comparisons=None, ca=False, outdir=None, outfile=None, subset=None, cluster_no=None, atoms=None, no_write=False, all_fuzzy=False, no_alt=False, find_ligands=True, transform_to_ref=True):
        """Creates template from conserved sites"""
        # Get reference as functional site (maybe make a separate method for retrieving reference from entry)
        if subset is None:
            subset = []
        try:
            reference = self._get_template_reference(ca=ca).get_functional_site(ca=ca)
        except TypeError:
            return
        # Calculate average coordinates
        avg = reference.copy(include_structure=True)
        for site in self.get_pdbsites(sane_only=True):
            if subset and site.id not in subset:
                continue
            if (not ca and not (site.is_conserved or site.is_conservative_mutation)) or \
                    site.id==reference.id or site.has_missing_functional_atoms:
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
        template = reference
        for site in self.get_pdbsites(sane_only=True):
            if subset and site.id not in subset:
                continue
            if (not ca and not (site.is_conserved or site.is_conservative_mutation)) or \
                    (site.id==reference.id or not ca and site.has_missing_functional_atoms):
                continue
            if 'nmr' in site.experimental_method.lower():
                continue
            rot, tran, rms, rms_all = avg.fit(site, transform=False, ca=ca)
            if rms_all < min_rms:
                min_rms = rms_all
                template = site
        # Fit template to reference or any other site with sane atoms
        if transform_to_ref:
            reference.fit(template, transform=True, ca=ca)
        # Bind subset to template
        template.subset = subset
        # Write template
        if no_write == False:
            self.write_template(template, comparisons, ca, template.subset, cluster_no, None, atoms, no_alt, all_fuzzy, outdir, outfile, '', find_ligands)
        return template

    def break_template(self, template, method='clusters', n_residues=3, max_distance=6):
        """Breaks template in smaller groups of arbitrary size (default=3).
        Returns a list of indeces of the residues of each group. Overlap is allowed."""
        if method=='combinations':
            groups = set()
            combis = list(combinations(range(template.size), n_residues))
            for combi in combis:
                residues = [template.residues[i] for i in combi]
                discard = False
                for i in residues:
                    for j in residues:
                        if i.get_distance(j, kind='com') > max_distance:
                            discard = True
                if not discard:
                    groups.add(combi)
            return list(groups)
        elif method=='clusters':
            # First calulate the centers of mass of each residue
            ids = {}
            coms = []
            for res in template:
                coms.append(res.structure.center_of_mass())
            coms = np.array(coms)
            # Calculate the number of groups to generate
            n_clusters = int(np.round(len(template.residues)/n_residues))
            # Do k-means clustering in template residues centers of mass
            kmeans = KMeans(n_clusters = n_clusters, random_state=0).fit(coms)
            # Generate groups around each cluster centre
            centers = kmeans.cluster_centers_
            groups = []
            for center in centers:
                dists = []
                for com in coms:
                    dists.append(np.linalg.norm(com-center))
                top = np.argsort(np.array(dists))[:n_residues]
                groups.append(top)
            return groups
        elif method=='functional':
            # Get free text annotation for each residue
            annotations = self._get_free_text_annotations()
            groups = defaultdict(list)
            for i, res in enumerate(template):
                for annotation in annotations[(res.reference_residue.chain, res.reference_residue.auth_resid)]:
                    groups[annotation].append(i)
            return groups

    def write_template(self, template, comparisons=None, ca=False, subset=None, cluster_no=None, residues=None, atoms=None, 
                       no_alt=False, all_fuzzy=False, outdir=None, outfile=None, annotation='', find_ligands=True):
        """
        Writes template coordinates and constraints in TESS/Jess format
        Args:
            outdir: Output directory
            outfile: Default 'mcsa_id.cluster_no.template.pdb' 
        """
        if type(subset) == list and len(subset)==0:
            subset = None
        if len(template.residues) < 3:
            return
        if not outdir:
            outdir = '.'
        if not outfile:
            cluster = 'all' if cluster_no is None else f'cluster_{cluster_no}'
            outfile = f'{outdir}/csa3d_{str(template.mcsa_id).zfill(4)}.{cluster}.{template.id}.template.pdb'
        size = len(subset) if subset else 'ALL'
        index = residues if residues else (i for i,j in enumerate(template.residues))
        index = '-'.join([str(i) for i in index])
        cluster_ligands, substrate_score, cofactor_score, ion_score, artefact_score = '','','','',''
        if find_ligands:
            cluster_ligands, substrate_score, cofactor_score, ion_score, artefact_score = self._get_cluster_ligands(template, subset, residues)
        with open(outfile, 'w') as o:
            remarks = (f'REMARK TEMPLATE\n'
                       f'REMARK CLUSTER {cluster_no}\n'
                       f'REMARK INDEX {index}\n'
                       f'REMARK REPRESENTING {size} CATALYTIC SITES\n'
                       f'REMARK ID {template.id}\n'
                       f'REMARK MCSA_ID {template.mcsa_id}\n'
                       f'REMARK PDB_ID {template.pdb_id}\n'
                       f'REMARK UNIPROT_ID {template.uniprot_id}\n'
                       f'REMARK EC {template.ec}\n'
                       f'REMARK ENZYME {template.enzyme}\n'
                       f'REMARK EXPERIMENTAL_METHOD {template.experimental_method}\n'
                       f'REMARK RESOLUTION {template.resolution}\n'
                       f'REMARK ORGANISM_NAME {template.organism_name}\n'
                       f'REMARK ORGANISM_ID {template.organism_id}\n'
                       f'REMARK ANNOTATION {annotation}\n'
                       f'REMARK LIGANDS_IN_CLUSTER {cluster_ligands}\n'
                       f'REMARK SUBSTRATE_SCORE {substrate_score}\n'
                       f'REMARK COFACTOR_SCORE {cofactor_score}\n'
                       f'REMARK ION_SCORE {ion_score}\n'
                       f'REMARK ARTEFACT_SCORE {artefact_score}')
            print(remarks, file=o)
            matchcodes, alt_residues = self._get_template_params(template, all_fuzzy)
            if no_alt:
                alt_residues = None
            dist_cutoffs = None
            if comparisons is not None:
                if not comparisons.empty:
                    dist_cutoffs = self.get_dist_cutoffs(comparisons, subset)
            for i, res in enumerate(template):
                if residues is not None and i not in residues:
                    continue
                if res.structure is not None:
                    resname = res.structure.get_resname()
                    if not res.is_standard:
                        resname = 'PTM'
                    if res.has_main_chain_function:
                        resname = 'ANY'
                    for atom in res.structure.get_atoms():
                        funcstring = '{}.{}'.format(resname if not ca else 'ANY', atom.get_id().upper())
                        if atoms is not None:
                            if atom.get_id().upper() not in atoms:
                                continue
                        else:
                            if funcstring not in RESIDUE_DEFINITIONS:
                                continue
                        matchcode = matchcodes.get(int(i), {}).get(atom.name, 0)
                        dist_cutoff = dist_cutoffs[i] if dist_cutoffs is not None else 0.0
                        alt = alt_residues[i] if alt_residues is not None else ''
                        pdb_line = '{:6}{:5d} {:<4}{}{:>3}{:>2}{:>4}{:>12.3f}{:>8.3f}{:>8.3f} {:<5.5} {:<5.2f}'.format(
                            'ATOM', matchcode, atom.name if len(atom.name) == 4 else ' {}'.format(atom.name), 'Z',
                            resname, res.structure.get_parent().get_id(), res.structure.get_id()[1],
                            atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2],
                            alt, dist_cutoff )
                        print(pdb_line, file=o)
            print('END', file=o)

    def _get_cluster_ligands(self, template, subset=None, residues=None):
        """Returns all ligands from the template and its homologues from the
        cluster it represents. Ligands must be within 5A away from all template residues"""
        ligands = []
        substrate_score = 0
        cofactor_score = 0
        ion_score = 0
        artefact_score = 0
        n = 0
        for site in self.pdbsites:
            has_substrate = 0
            has_cofactor = 0
            has_ion = 0
            has_artefact = 0
            if subset is not None and site.id not in subset:
                continue
            n+=1
            for ligand in site.ligands:
                reject = False
                for i, res in enumerate(site.residues):
                    if residues is not None and i not in residues:
                        continue
                    d = round(res.get_distance(ligand, kind='min'), 2)
                    if d>6 or np.isnan(d):
                        reject = True
                if not reject:
                    ligands.append(f'{ligand.resname}:{ligand.type}:{d}')
                    if ligand.type == 'Substrate (polymer)':
                        has_substrate = 1
                    if ligand.type == 'Substrate (non-polymer)' and (ligand.similarity is not None and ligand.similarity>=0.5):
                        has_substrate = 1
                    if ligand.type == 'Co-factor (non-ion)':
                        has_cofactor = 1
                    if ligand.type in ('Co-factor (ion)', 'Ion'):
                        has_ion = 1
                    if ligand.type == 'Artefact':
                        has_artefact = 1
            substrate_score += has_substrate
            cofactor_score += has_cofactor
            ion_score += has_ion
            artefact_score += has_artefact
        if n==0:
            n=1
        return ';'.join(ligands), np.round(substrate_score/n, 2), \
               np.round(cofactor_score/n, 2), np.round(ion_score/n, 2), np.round(artefact_score/n, 2)

    def _get_template_params(self, template, all_fuzzy=False):
        alt_residues, backbone_fuzzy = self._get_alt_residues(template, all_fuzzy)
        matchnumbers = self._get_matchnumbers(template, backbone_fuzzy)
        return matchnumbers, alt_residues
        
    def _get_alt_residues(self, template, all_fuzzy=False):
        """Returns a string of alternative residues to be matched in the template"""
        observed = defaultdict(set)
        expected = defaultdict(set)
        intersection = defaultdict(set)
        backbone_fuzzy = {}
        # Get alternative residues in each position from sequence alignment
        for site in template.parent_entry.pdbsites:
            if site.is_conserved or site.is_conservative_mutation:
                for i, res in enumerate(site):
                    observed[i].add(AA_3TO1[res.resname])
        # Compare expected with observed substitutions
        ref = template.reference_site
        for i, res in enumerate(ref):
            expected[i] = set(list(TEMPLATE_FUZZY_RESIDUES.get(res.resname.upper(), '')))
            expected[i].add(AA_3TO1[res.resname])
        for i in expected:
            backbone_fuzzy[i] = True
            if all_fuzzy:
                intersection[i] = ''.join(list(expected[i])).replace(' ','')
            else:
                intersection[i] = ''.join(list(expected[i].intersection(observed[i]))).replace(' ','')
                res = site.residues[i]
                if res.has_main_chain_function or 'ptm' in res.funclocs or not res.is_standard:
                    if len(observed[i])<=5:
                        intersection[i] = ''.join(list(observed[i])).replace(' ','')
                        backbone_fuzzy[i] = False
        return intersection, backbone_fuzzy

    def _get_matchnumbers(self, template, backbone_fuzzy, ca=False):
        """Returns the Jess match option number for each template atom"""
        matchnumbers = defaultdict(dict)
        for i, res in enumerate(template):
            resname = res.resname.upper()
            if 'ptm' in res.funclocs or not res.is_standard:
                resname = 'PTM'
            for funcloc in res.funclocs:
                if 'main' in funcloc:
                    resname = 'ANY'
            if not res.structure:
                continue
            for atom in res.structure:
                if not backbone_fuzzy[i]:
                    matchnumbers[i][atom.name] = 0
                    continue
                if ca:
                    if atom.name in TEMPLATE_FUZZY_ATOMS['ANY']:
                        if resname in ['ANY', 'PTM']:
                            matchnumbers[i][atom.name] = 100
                        else:
                            matchnumbers[i][atom.name] = 0
                else:
                    if atom.name in TEMPLATE_FUZZY_ATOMS[resname]:
                        matchnumbers[i][atom.name] = TEMPLATE_FUZZY_ATOMS[resname][atom.name][0]
        return matchnumbers

    def get_dist_cutoffs(self, comparisons, subset=None):
        """Returns the distance threshold weight for each residue"""
        if type(comparisons) == str:
            comparisons = pd.read_csv(open(comparisons, 'r'))
        if subset:
            if len(subset) <= 1:
                return
            comparisons = comparisons.query('p_id in @subset and q_id in @subset')
        per_res = comparisons.loc[~pd.isnull(comparisons['per_res_rms']), 'per_res_rms'].astype('str').str.split(';', expand=True).astype('float32')
        try:
            cutoffs = dict(per_res.describe().loc[['mean', 'std']].sum())
        except ValueError:
            return
        return cutoffs

    def _get_free_text_annotations(self):
        """Gets free text role annotations for each reference residue
        from the M-CSA API file and returns them as a list, preserving
        the order of the residues as chiral id in the reference catalytic site"""
        info = json.load(open(REF_RES_INFO))
        annotations = defaultdict(list)
        for res in info:
            if res['mcsa_id'] == self.mcsa_id:
                chain, resid = res['residue_chains'][0]['assembly_chain_name'], res['residue_chains'][0]['auth_resid']
                annotations[(chain, resid)].append(res['main_annotation'])
        return annotations

