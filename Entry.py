"""M-CSA entry representation and analysis.

Provides the Entry class for manipulating enzyme catalytic sites from both
PDB structures and UniProt sequences. Includes RMSD calculation, clustering,
and template generation methods.
"""

import logging
import os
import json
import numpy as np
import pandas as pd
from copy import deepcopy
from collections import defaultdict
from itertools import combinations
from typing import Dict, List, Optional, Set, Tuple, Generator, Any
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import linkage, cut_tree
from sklearn.cluster import KMeans
from ._vendor.dendrogram_cut import DendrogramCut
from .config import UNI2PDB, REF_RES_INFO
from .residue_definitions import (
    RESIDUE_DEFINITIONS, TEMPLATE_FUZZY_RESIDUES, TEMPLATE_FUZZY_ATOMS, AA_3TO1
)
from .PdbSite import PdbSite
from .UniSite import UniSite

logger = logging.getLogger(__name__)


class Entry:
    """M-CSA entry containing catalytic sites from PDB and UniProt.

    Represents a single M-CSA catalytic mechanism entry, containing both
    PDB structure sites and UniProt sequence sites. Provides analysis tools
    including RMSD clustering, template creation, and site comparison.

    Attributes:
        mcsa_id: Unique M-CSA identifier
        info: Metadata dictionary (EC, enzyme name, organism, etc.)
        pdbsites: List of PdbSite objects (3D structures)
        unisites: List of UniSite objects (sequences)
        pdb_ids: Set of unique PDB IDs in this entry
        pdbsites_dict: Dictionary mapping site ID to PdbSite
        unisites_dict: Dictionary mapping site ID to UniSite

    Example:
        >>> entry = Entry(mcsa_id=1)
        >>> len(entry.pdbsites)
        0
    """

    def __init__(self, mcsa_id: Optional[int] = None) -> None:
        """Initialize an M-CSA entry.

        Args:
            mcsa_id: M-CSA entry identifier. Defaults to None.
        """
        self.mcsa_id: Optional[int] = mcsa_id
        self.info: Optional[Dict[str, Any]] = None
        self.unisites: List[UniSite] = []
        self.unisites_dict: Dict[str, UniSite] = {}
        self.pdbsites: List[PdbSite] = []
        self.pdbsites_dict: Dict[str, PdbSite] = {}
        self.pdb_ids: Set[str] = set()

    @property
    def has_multiple_pdb_refs(self) -> bool:
        """Check if entry has more than one PDB reference structure."""
        refs = set()
        for site in self.pdbsites:
            if site.is_reference:
                refs.add(site.pdb_id)
        return len(refs) > 1

    @property
    def has_multiple_uniprot_refs(self) -> bool:
        """Check if entry has more than one UniProt reference sequence."""
        refs = set()
        for site in self.unisites:
            if site.is_reference:
                refs.add(site.uniprot_id)
        return len(refs) > 1

    def add(self, site: Any) -> bool:
        """Add a catalytic site (PdbSite or UniSite) to this entry.

        Automatically routes site to appropriate list and updates internal
        mappings between PDB and UniProt representations.

        Args:
            site: PdbSite or UniSite object to add

        Returns:
            True if successful, False if wrong type

        Raises:
            TypeError: If site is not PdbSite or UniSite (future versions)
        """
        if isinstance(site, UniSite):
            self.unisites.append(site)
            self.unisites_dict[site.id] = site
        elif isinstance(site, PdbSite):
            self.pdbsites.append(site)
            self.pdbsites_dict[site.id] = site
            self.pdb_ids.add(site.pdb_id)
        else:
            logger.error(f"Cannot add {type(site).__name__}. Expected PdbSite or UniSite.")
            return False

        self.update_mapped_sites()
        return True

    def get_unisite(self, _id: str) -> Optional[UniSite]:
        """Get UniSite by ID.

        Args:
            _id: UniProt accession ID

        Returns:
            UniSite object or None if not found
        """
        return self.unisites_dict.get(_id)

    def get_pdbsite(self, _id: str) -> Optional[PdbSite]:
        """Get PdbSite by ID.

        Args:
            _id: Site identifier

        Returns:
            PdbSite object or None if not found
        """
        return self.pdbsites_dict.get(_id)

    def get_unisites(self) -> Generator[UniSite, None, None]:
        """Get all UniSite objects.

        Yields:
            UniSite: Each UniSite in this entry
        """
        yield from self.unisites

    def get_pdbsites(
        self, pdb_id: Optional[str] = None, sane_only: bool = False
    ) -> Generator[PdbSite, None, None]:
        """Get PdbSite objects, optionally filtered.

        Args:
            pdb_id: Filter by PDB ID. Defaults to None (all sites).
            sane_only: Only return sites marked as 'sane'. Defaults to False.

        Yields:
            PdbSite: Matching sites
        """
        if pdb_id is None:
            for site in self.pdbsites:
                if sane_only and not site.is_sane:
                    continue
                yield site
        else:
            result = [
                val for key, val in self.pdbsites_dict.items() if pdb_id in key
            ]
            for site in result:
                if sane_only and not site.is_sane:
                    continue
                yield site

    def update_mapped_sites(self) -> bool:
        """Update mappings between PDB and UniProt sites.

        Uses SIFTS mapping to link UniSites to their corresponding PdbSites,
        enabling cross-representation analysis.

        Returns:
            True
        """
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

    def all_vs_all(
        self, sane_only: bool = False
    ) -> Generator[Tuple[PdbSite, PdbSite], None, None]:
        """Generate all pairwise combinations of PDB sites.

        Args:
            sane_only: Only include sites marked as 'sane'. Defaults to False.

        Yields:
            Tuple of two PdbSite objects (p, q) where p != q
        """
        seen: Set[Tuple[str, str]] = set()
        for p in self.pdbsites:
            if sane_only and not p.is_sane:
                continue
            for q in self.pdbsites:
                if sane_only and not q.is_sane:
                    continue
                if p == q or (q.id, p.id) in seen:
                    continue
                seen.add((p.id, q.id))
                yield p, q

    def rmsd_matrix(
        self, sitelist: List[PdbSite], ca: bool = False
    ) -> pd.DataFrame:
        """Calculate pairwise RMSD matrix for catalytic sites.

        Computes structural similarity between all pairs of sites using
        spatial RMSD (root-mean-square-deviation).

        Args:
            sitelist: List of PdbSite objects to compare
            ca: Use CA (alpha carbon) only. Defaults to False (all atoms).

        Returns:
            DataFrame with RMSD values (symmetric, NaN on diagonal)

        Raises:
            ValueError: If insufficient sites or missing coordinates
        """
        rmsds: List[float] = []
        ids: List[str] = []
        seen: Set[Tuple[str, str]] = set()

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
                except Exception as e:
                    logger.debug(f"Weighted fit failed for {p.id}-{q.id}: {e}, trying unweighted")
                    _, _, _, rms_all = p.fit(q, weighted=False, ca=ca, cycles=10, cutoff=5)
                rmsds.append(rms_all)

        rmsds_array = np.array(rmsds, dtype='float32')
        matrix = pd.DataFrame(squareform(rmsds_array), columns=ids, index=ids)
        return matrix

    def clustering(
        self,
        matrix: pd.DataFrame,
        height: Optional[float] = None,
        get_linkage_matrix: bool = False,
    ) -> Union[Dict[int, List[str]], Tuple[Dict[int, List[str]], np.ndarray]]:
        """Perform hierarchical clustering on RMSD matrix.

        Groups catalytic sites based on structural similarity using average
        linkage hierarchical clustering.

        Args:
            matrix: RMSD distance matrix (from rmsd_matrix())
            height: Cutoff height for tree cutting. If None, uses automatic
                (70% of max distance). Defaults to None.
            get_linkage_matrix: Also return the linkage matrix. Defaults to False.

        Returns:
            Dictionary mapping cluster ID to list of site IDs. Optionally includes
            the linkage matrix as second element of tuple.

        Example:
            >>> entry = Mcsa().entries[1]
            >>> sites = list(entry.get_pdbsites(sane_only=True))
            >>> matrix = entry.rmsd_matrix(sites)
            >>> clusters = entry.clustering(matrix)
            >>> len(clusters)
            3
        """
        ids = list(matrix.index)
        if len(ids) < 2:
            return {0: ids}

        # Compute linkage matrix
        Z = linkage(squareform(matrix), method='average')

        # Cut tree at automatic height
        if height is None:
            height = 0.7 * max(Z[:, 2])

        # Get clusters index
        clusters = cut_tree(Z, height=height).flatten()
        # Group sites to clusters
        cluster_dict: Dict[int, List[str]] = defaultdict(list)
        for cluster in np.unique(clusters):
            for i, site_id in zip(clusters, ids):
                if i == cluster:
                    cluster_dict[cluster].append(self.get_pdbsite(site_id).id)
        if get_linkage_matrix:
            return cluster_dict, Z
        return cluster_dict

    def clustering_bayesian(
        self,
        matrix: pd.DataFrame,
        k_max: int = 5,
        l: float = 1.0,
        plot_outfile: Optional[str] = None,
    ) -> Dict[int, List[str]]:
        """Perform hierarchical clustering with Bayesian automatic cutting.

        Uses PAC-Bayesian statistics to automatically determine optimal
        number of clusters. More sophisticated than fixed-height cutting.

        Args:
            matrix: RMSD distance matrix (from rmsd_matrix())
            k_max: Maximum number of clusters to test. Defaults to 5.
            l: Lambda parameter for Bayesian cut. Defaults to 1.0.
            plot_outfile: Optional path to save dendrogram plot. Defaults to None.

        Returns:
            Dictionary mapping cluster ID to list of site IDs

        Raises:
            ValueError: If matrix is invalid
        """
        cluster_dict: Dict[int, List[str]] = defaultdict(list)
        ids = list(matrix.index)
        matrix_array = np.array(matrix)
        try:
            model = DendrogramCut(k_max=k_max, method='average').fit(matrix_array)
        except ValueError as e:
            logger.warning(f"Bayesian clustering failed: {e}. Using single cluster.")
            return {0: ids}
        try:
            k = model.pac_bayesian_cut(lambda_=l)
        except RecursionError:
            logger.warning("Recursion error in Bayesian cut. Using single cluster.")
            return {0: ids}
        clusters = model.get_cluster_label(k=k)
        if plot_outfile:
            fig = model.heatmap_with_dendrogram_plot(k=k)
            fig.write_image(plot_outfile)
            logger.info(f"Saved dendrogram plot to {plot_outfile}")
        for i, cluster_no in enumerate(clusters):
            cluster_dict[cluster_no].append(ids[i])
        return cluster_dict

    def _get_template_reference(self, ca: bool = False) -> Optional[PdbSite]:
        """Get a reference site for template generation.

        Selects the best reference site for fitting other sites into,
        prioritizing reference structures over conserved homologues.

        Args:
            ca: Use CA (alpha carbon) only. Defaults to False.

        Returns:
            Best available PdbSite for template reference, or None if none found
        """
        ref: Optional[PdbSite] = None
        for site in self.get_pdbsites():
            if site.is_reference:
                if ca or (not ca and not site.has_missing_functional_atoms):
                    ref = site
                    break
        if not ref:
            for site in self.get_pdbsites(sane_only=True):
                if site.is_conserved or site.is_conservative_mutation:
                    if ca or (not ca and not site.has_missing_functional_atoms):
                        ref = site
                        break
        return ref

    def create_template(
        self,
        comparisons: Optional[Union[str, pd.DataFrame]] = None,
        ca: bool = False,
        outdir: Optional[str] = None,
        outfile: Optional[str] = None,
        subset: Optional[List[str]] = None,
        cluster_no: Optional[int] = None,
        atoms: Optional[List[str]] = None,
        no_write: bool = False,
        all_fuzzy: bool = False,
        no_alt: bool = False,
    ) -> Optional[PdbSite]:
        """Create a template structure from conserved catalytic sites.

        Generates an averaged template structure representing the consensus
        active site geometry across conserved homologues.

        Args:
            comparisons: Comparison data for per-residue constraints (optional)
            ca: Use CA atoms only. Defaults to False.
            outdir: Output directory for template file. Defaults to current dir.
            outfile: Output filename. Defaults to auto-generated name.
            subset: Limit to specific site IDs. Defaults to None (all conserved).
            cluster_no: Cluster identifier for output filename. Defaults to None.
            atoms: Specific atoms to include. Defaults to None (functional atoms).
            no_write: Don't write output file. Defaults to False.
            all_fuzzy: Allow all fuzzy residue matches. Defaults to False.
            no_alt: Exclude alternative residues. Defaults to False.

        Returns:
            PdbSite object representing the template, or None if failed

        Example:
            >>> entry = Mcsa().entries[1]
            >>> template = entry.create_template(outdir="./templates")
            >>> template.pdb_id
            '1abc'
        """
        if subset is None:
            subset = []
        try:
            reference = self._get_template_reference(ca=ca).get_functional_site(ca=ca)
        except (TypeError, AttributeError) as e:
            logger.error(f"Could not get template reference: {e}")
            return None

        # Calculate average coordinates
        avg = reference.copy(include_structure=True)
        for site in self.get_pdbsites(sane_only=True):
            if subset and site.id not in subset:
                continue
            if (not ca and not (site.is_conserved or site.is_conservative_mutation)) or \
                    site.id == reference.id or site.has_missing_functional_atoms:
                continue
            funcsite = site.get_functional_site(ca=ca)
            avg.fit(funcsite, transform=True, ca=ca)
            site_coords = np.array([atom.get_coord() for res in funcsite for atom in res.structure])
            avg_coords = np.array([atom.get_coord() for res in avg for atom in res.structure])
            avg_coords = np.mean([avg_coords, site_coords], axis=0)
            for atom, coord in zip([atom for res in avg for atom in res.structure], avg_coords):
                atom.set_coord(coord)

        # Find representative site
        min_rms = 999.0
        template = reference
        for site in self.get_pdbsites(sane_only=True):
            if subset and site.id not in subset:
                continue
            if (not ca and not (site.is_conserved or site.is_conservative_mutation)) or \
                    (site.id == reference.id or (not ca and site.has_missing_functional_atoms)):
                continue
            rot, tran, rms, rms_all = avg.fit(site, transform=False, ca=ca)
            if rms_all < min_rms:
                min_rms = rms_all
                template = site

        # Fit template to reference or any other site with sane atoms
        reference.fit(template, transform=True, ca=ca)
        template.subset = subset

        # Write template if requested
        if not no_write:
            self.write_template(
                template, comparisons, ca, template.subset, cluster_no,
                None, atoms, no_alt, all_fuzzy, outdir, outfile
            )
        return template

    def break_template(
        self,
        template: PdbSite,
        method: str = 'clusters',
        n_residues: int = 3,
        max_distance: float = 6.0,
    ) -> Union[List[List[int]], Dict[str, List[int]]]:
        """Break template into smaller residue groups.

        Divides a large template into smaller, manageable subsets for
        more focused analysis or docking.

        Args:
            template: Template PdbSite to break
            method: Grouping method - 'clusters', 'combinations', or 'functional'.
                Defaults to 'clusters'.
            n_residues: Target residues per group. Defaults to 3.
            max_distance: Maximum distance between residues in group (Å).
                Defaults to 6.0.

        Returns:
            List of residue index groups, or dict if using 'functional' method

        Example:
            >>> template = entry.create_template()
            >>> groups = entry.break_template(template, method='clusters', n_residues=3)
        """
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

    def write_template(
        self,
        template: PdbSite,
        comparisons: Optional[pd.DataFrame] = None,
        ca: bool = False,
        subset: Optional[List[str]] = None,
        cluster_no: Optional[int] = None,
        residues: Optional[List[int]] = None,
        atoms: Optional[List[str]] = None,
        no_alt: bool = False,
        all_fuzzy: bool = False,
        outdir: Optional[str] = None,
        outfile: Optional[str] = None,
        annotation: str = '',
    ) -> None:
        """Write template in PDB/TESS format with constraints.

        Exports template structure with geometric constraints for use in
        docking software (TESS/Jess).

        Args:
            template: Template PdbSite to write
            comparisons: RMSD comparison data for constraints. Defaults to None.
            ca: Use CA atoms only. Defaults to False.
            subset: Limit to specific site IDs. Defaults to None.
            cluster_no: Cluster identifier. Defaults to None.
            residues: Specific residue indices to include. Defaults to None.
            atoms: Specific atoms to include. Defaults to None.
            no_alt: Exclude alternative residue options. Defaults to False.
            all_fuzzy: Allow all fuzzy matches. Defaults to False.
            outdir: Output directory. Defaults to current directory.
            outfile: Output filename. Auto-generated if None.
            annotation: Additional annotation comment. Defaults to empty.

        Returns:
            None (writes file to disk)
        """
        if isinstance(subset, list) and len(subset) == 0:
            subset = None
        if len(template.residues) < 3:
            logger.warning("Template has fewer than 3 residues, skipping write")
            return
        if not outdir:
            outdir = '.'
        if not outfile:
            cluster = 'all' if cluster_no is None else f'cluster_{cluster_no}'
            outfile = f'{outdir}/csa3d_{str(template.mcsa_id).zfill(4)}.{cluster}.{template.id}.template.pdb'

        size = len(subset) if subset else 'ALL'
        cluster_ligands, substrate_score, cofactor_score, ion_score, artefact_score = self._get_cluster_ligands(
            template, subset, residues
        )
        with open(outfile, 'w') as o:
            remarks = (f'REMARK TEMPLATE\n'
                       f'REMARK CLUSTER {cluster_no}\n'
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

    def _get_cluster_ligands(
        self,
        template: PdbSite,
        subset: Optional[List[str]] = None,
        residues: Optional[List[int]] = None,
    ) -> Tuple[str, float, float, float, float]:
        """Get ligands associated with template cluster.

        Returns ligands from the template and its homologues that are within
        5Å of any catalytic residue.

        Args:
            template: Template PdbSite
            subset: Limit to specific site IDs. Defaults to None.
            residues: Limit to specific residue indices. Defaults to None.

        Returns:
            Tuple of (ligands_str, substrate_score, cofactor_score,
                     ion_score, artefact_score)
        """
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

    def _get_template_params(
        self, template: PdbSite, all_fuzzy: bool = False
    ) -> Tuple[Dict[int, Dict[str, int]], Dict[int, str]]:
        """Get template parameters (matchcodes and alternative residues).

        Args:
            template: Template PdbSite
            all_fuzzy: Allow all fuzzy matches. Defaults to False.

        Returns:
            Tuple of (matchcodes dict, alt_residues dict)
        """
        alt_residues, backbone_fuzzy = self._get_alt_residues(template, all_fuzzy)
        matchnumbers = self._get_matchnumbers(template, backbone_fuzzy)
        return matchnumbers, alt_residues

    def _get_alt_residues(
        self, template: PdbSite, all_fuzzy: bool = False
    ) -> Tuple[Dict[int, str], Dict[int, bool]]:
        """Get alternative residues for fuzzy matching.

        Determines which amino acids can be matched at each position based on
        observed substitutions in homologous sites.

        Args:
            template: Template PdbSite
            all_fuzzy: Allow all fuzzy matches. Defaults to False.

        Returns:
            Tuple of (alternatives dict, backbone_fuzzy dict)
        """
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

    def _get_matchnumbers(
        self,
        template: PdbSite,
        backbone_fuzzy: Dict[int, bool],
        ca: bool = False,
    ) -> Dict[int, Dict[str, int]]:
        """Get Jess/TESS match codes for template atoms.

        Determines the matching flexibility for each atom in the template.

        Args:
            template: Template PdbSite
            backbone_fuzzy: Dict indicating backbone fuzzy matching per position
            ca: Use CA atoms only. Defaults to False.

        Returns:
            Dictionary mapping residue index to atom name to match code
        """
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

    def get_dist_cutoffs(
        self,
        comparisons: Union[str, pd.DataFrame],
        subset: Optional[List[str]] = None,
    ) -> Optional[Dict[int, float]]:
        """Get distance cutoffs per residue from comparison data.

        Calculates per-residue distance constraints from RMSD comparisons.

        Args:
            comparisons: Path to CSV or DataFrame with comparison results
            subset: Limit to specific site IDs. Defaults to None.

        Returns:
            Dictionary mapping residue index to distance cutoff, or None if failed
        """
        if isinstance(comparisons, str):
            comparisons = pd.read_csv(open(comparisons, 'r'))
        if subset:
            if len(subset) <= 1:
                return None
            comparisons = comparisons.query('p_id in @subset and q_id in @subset')
        per_res = comparisons.loc[
            ~pd.isnull(comparisons['per_res_rms']), 'per_res_rms'
        ].astype('str').str.split(';', expand=True).astype('float32')
        try:
            cutoffs = dict(per_res.describe().loc[['mean', 'std']].sum())
            return cutoffs
        except ValueError as e:
            logger.warning(f"Could not calculate distance cutoffs: {e}")
            return None

    def _get_free_text_annotations(
        self,
    ) -> Dict[Tuple[str, int], List[str]]:
        """Get functional annotations for reference residues.

        Retrieves free-text role annotations from M-CSA API data for each
        reference catalytic residue.

        Returns:
            Dictionary mapping (chain, residue_id) to list of annotations
        """
        try:
            info = json.load(open(REF_RES_INFO))
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.error(f"Could not load reference residue info: {e}")
            return {}

        annotations: Dict[Tuple[str, int], List[str]] = defaultdict(list)
        for res in info:
            if res['mcsa_id'] == self.mcsa_id:
                chain = res['residue_chains'][0]['assembly_chain_name']
                resid = res['residue_chains'][0]['auth_resid']
                annotations[(chain, resid)].append(res['main_annotation'])
        return annotations


__all__ = ['Entry']