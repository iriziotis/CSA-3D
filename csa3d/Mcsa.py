"""M-CSA database access and entry building.

Provides the main Mcsa class for loading and building enzyme catalytic site
entries from the Mechanism and Catalytic Site Atlas (M-CSA) database.
"""

import logging
import os
import json
from collections import defaultdict
from glob import glob
from typing import Dict, List, Optional, Union
from natsort import natsorted

from . import datasets
from .data_setup import fetch_assembly
from .Entry import Entry
from .PdbSite import PdbSite
from .UniSite import UniSite
from .PdbResidue import PdbResidue
from .UniResidue import UniResidue

logger = logging.getLogger(__name__)


class Mcsa:
    """Top-level database class for M-CSA entries.

    Contains Entry objects and provides methods to load and build entries
    from M-CSA JSON data. Manages catalytic sites from PDB structures and
    UniProt sequences.

    Attributes:
        entries: Dictionary mapping M-CSA IDs to Entry objects
        mcsa_entry_info: Metadata for each M-CSA entry
        catalytic_residue_info: Catalytic residue definitions from M-CSA
        pdb_residues: Nested dict of residues organized by entry and PDB ID
        uni_residues: Nested dict of UniProt residues by entry and accession
        ref_pdb_residues: Reference PDB residues (highest quality structures)
        ref_uni_residues: Reference UniProt residues (canonical sequences)

    Example:
        >>> db = Mcsa()
        >>> db.build(entry_ids=1)
        >>> entry = db.entries[1]
        >>> len(entry.pdbsites)
        42
    """

    def __init__(self) -> None:
        """Initialize the M-CSA database."""
        self.entries: Dict[int, Entry] = {}
        self.mcsa_entry_info: Dict[int, dict] = self._get_mcsa_entry_info(datasets.MCSA_ENTRY_INFO)
        self.catalytic_residue_info: Dict[int, list] = self._get_residue_info(datasets.CAT_RES_INFO)
        self.pdb_assemblies: Dict[str, int] = self._get_pdb_assemblies(self.catalytic_residue_info)
        self.pdb_residues: Dict[int, Dict[str, List]] = defaultdict(lambda: defaultdict(list))
        self.uni_residues: Dict[int, Dict[str, List]] = defaultdict(lambda: defaultdict(list))
        self.ref_pdb_residues: Dict[int, List] = defaultdict(list)
        self.ref_uni_residues: Dict[int, List] = defaultdict(list)

    def __iter__(self):
        """Iterate over all Entry objects in the database."""
        yield from self.entries.values()

    def build(
        self,
        entry_ids: Optional[Union[int, List[int]]] = None,
        annotate: bool = True,
        redundancy_cutoff: Optional[float] = None,
        resolution_cutoff: Optional[float] = None,
        verbose: bool = False,
        no_sites: bool = False,
    ) -> bool:
        """Build M-CSA entries with catalytic sites.

        Constructs Entry objects containing PDB and UniProt catalytic sites
        from M-CSA JSON data.

        Args:
            entry_ids: M-CSA ID(s) to build. Can be int, list of ints, or None
                (build all). Defaults to None.
            annotate: Extract metadata (EC, organism, resolution, etc.).
                Defaults to True.
            redundancy_cutoff: Drop sites within this RMSD (Å) of another site in
                the same structure. Defaults to None (only identical sites are dropped).
            resolution_cutoff: Skip non-reference structures with resolution worse
                than this (Å), or with no resolution. Defaults to None (no filter).
            verbose: Print progress for each site. Defaults to False.
            no_sites: Skip site building (only build residues). Defaults to False.

        Returns:
            True if successful, False if entry not found or not yet implemented.

        Raises:
            KeyError: If entry_id not in M-CSA database

        Example:
            >>> db = Mcsa()
            >>> db.build(entry_ids=[1, 2, 3], verbose=True)
            True
        """
        if entry_ids is None:
            entry_ids = list(self.catalytic_residue_info.keys())
        elif not isinstance(entry_ids, list):
            entry_ids = [entry_ids]

        for entry_id in entry_ids:
            if entry_id not in self.catalytic_residue_info:
                logger.warning(f"Entry {entry_id} not found in M-CSA database")
                return False

            self.entries[entry_id] = Entry(entry_id)
            self.entries[entry_id].info = self.mcsa_entry_info[entry_id]
            self._build_pdb_residues(entry_id, resolution_cutoff)
            self._build_uniprot_residues(entry_id)

            # Check for unimplemented cases
            pdb_refs = set(r.pdb_id for r in self.ref_pdb_residues[entry_id])
            if len(pdb_refs) > 1:
                logger.error(
                    f"Entry {entry_id} has {len(pdb_refs)} PDB reference structures. "
                    "Multiple references not yet implemented."
                )
                return False

            uni_refs = set(r.uniprot_id for r in self.ref_uni_residues[entry_id])
            if len(uni_refs) > 1:
                logger.error(
                    f"Entry {entry_id} has {len(uni_refs)} UniProt reference sequences. "
                    "Multiple references not yet implemented."
                )
                return False

            if not no_sites:
                self._build_pdb_sites(entry_id, annotate, redundancy_cutoff, verbose)
                self._build_uniprot_sites(entry_id)

        return True

    def add(self, entry: Entry) -> bool:
        """Add an Entry to the database.

        Args:
            entry: Entry object to add

        Returns:
            True if successful, False if not an Entry

        Raises:
            TypeError: If entry is not an Entry instance
        """
        if not isinstance(entry, Entry):
            logger.error(f"Cannot add {type(entry).__name__}. Expected Entry object.")
            return False
        self.entries[entry.mcsa_id] = entry
        logger.debug(f"Added entry {entry.mcsa_id} to database")
        return True

    def update_pdbs(self, entry_ids: Optional[Union[int, List[int]]] = None) -> None:
        """Pre-download assembly mmCIF files from PDBe.

        ``build()`` already downloads missing assemblies on demand; use this to
        fetch everything up front (e.g. before working offline). Files already
        present in the data directory are skipped.

        Args:
            entry_ids: M-CSA ID(s) whose structures to fetch. Defaults to all.
        """
        if entry_ids is None:
            entry_ids = list(self.catalytic_residue_info)
        elif not isinstance(entry_ids, list):
            entry_ids = [entry_ids]
        pdb_ids = {
            res_chain['pdb_id']
            for entry_id in entry_ids
            for res in self.catalytic_residue_info.get(entry_id, [])
            for res_chain in res['residue_chains']
        }
        for pdb_id in sorted(pdb_ids):
            self._get_cif_path(pdb_id)

    # Properties

    @property
    def nofentries(self) -> int:
        """Number of entries in the database."""
        return len(self.entries)

    # Private methods

    def _build_pdb_residues(
        self, entry_id: int, resolution_cutoff: Optional[float] = None
    ) -> None:
        """Build PdbResidue objects from M-CSA JSON data.

        Constructs PdbResidue instances from raw catalytic residue information.
        Extracts all equivalent residues from identical chains in biological
        assemblies and applies resolution filtering if specified.

        Args:
            entry_id: M-CSA entry ID
            resolution_cutoff: Only include structures with resolution ≤ this value (Å).
        """
        reference_residue = None
        for chiral_id, json_res in enumerate(self.catalytic_residue_info[entry_id]):
            for residue in PdbResidue.from_json(json_res, chiral_id):
                # Resolution cutoff
                if resolution_cutoff and not residue.is_reference:
                    try:
                        resolution = float(datasets.PDBID_RESOLUTION.get(residue.pdb_id, -1.0))
                    except ValueError:
                        resolution = -1.0
                    if resolution == -1.0 or resolution > resolution_cutoff:
                        continue
                if residue.is_reference:
                    reference_residue = residue
                    self.ref_pdb_residues[residue.mcsa_id].append(residue)
                else:
                    self.pdb_residues[residue.mcsa_id][residue.pdb_id].append(residue)
                residue.reference_residue = reference_residue

    def _build_uniprot_residues(self, entry_id: int) -> None:
        """Build UniResidue objects from M-CSA JSON data.

        Constructs UniResidue instances representing catalytic residues in
        UniProt sequences, including homologous residues and their alignment
        to the reference sequence.

        Args:
            entry_id: M-CSA entry ID
        """
        reference_residue = None
        for chiral_id, json_res in enumerate(self.catalytic_residue_info[entry_id]):
            for residue in UniResidue.from_json(json_res, chiral_id):
                if residue.is_reference:
                    reference_residue = residue
                    self.ref_uni_residues[residue.mcsa_id].append(residue)
                else:
                    self.uni_residues[residue.mcsa_id][residue.uniprot_id].append(residue)
                residue.reference_residue = reference_residue

    def _build_pdb_sites(
        self,
        entry_id: int,
        annotate: bool = True,
        redundancy_cutoff: Optional[float] = None,
        verbose: bool = False,
    ) -> None:
        """Build PdbSite objects from PdbResidue lists.

        Constructs PdbSite instances using geometric distance criteria.
        Automatically detects homologous sites in different structures
        and groups them by catalytic equivalence.

        Args:
            entry_id: M-CSA entry ID
            annotate: Extract structural annotations (resolution, method, etc.).
                Defaults to True.
            redundancy_cutoff: Drop sites within this RMSD (Å) of another site in
                the same structure. Defaults to None (only identical sites are dropped).
            verbose: Print each site as it's built. Defaults to False.
        """
        try:
            ref_pdb_id = self.ref_pdb_residues[entry_id][0].pdb_id
        except (IndexError, KeyError):
            logger.warning(f"Entry {entry_id}: No reference PDB structure found")
            return

        reference_site = PdbSite.build_reference(
            self.ref_pdb_residues[entry_id],
            self.entries[entry_id],
            self._get_cif_path(ref_pdb_id),
            annotate,
        )
        self.entries[entry_id].add(reference_site)

        for pdb_id, pdb_residues in self.pdb_residues[entry_id].items():
            for site in PdbSite.build_all(
                pdb_residues,
                reference_site,
                self.entries[entry_id],
                self._get_cif_path(pdb_id),
                annotate,
                redundancy_cutoff,
            ):
                if verbose:
                    logger.info(f"Built site: {site.id}")
                self.entries[entry_id].add(site)

    def _build_uniprot_sites(self, entry_id: int) -> None:
        """Build UniSite objects from UniResidue lists.

        Constructs UniSite instances representing catalytic sites at the
        sequence level across different UniProt entries.

        Args:
            entry_id: M-CSA entry ID
        """
        reference_site = UniSite.build(self.ref_uni_residues[entry_id])
        reference_site.parent_entry = self.entries[entry_id]
        self.entries[entry_id].add(reference_site)
        for uniprot_id, uniprot in self.uni_residues[entry_id].items():
            site = UniSite.build(uniprot, reference_site)
            site.parent_entry = self.entries[entry_id]
            self.entries[entry_id].add(site)

    # Static methods

    @staticmethod
    def _get_mcsa_entry_info(jsons_wildcard: str) -> Dict[int, dict]:
        """Load M-CSA entry metadata from JSON files.

        Loads basic information for all M-CSA entries from JSON files matching
        the provided wildcard pattern. Entries are indexed by M-CSA ID.

        Args:
            jsons_wildcard: Glob pattern for JSON files (e.g., "data/*.json")

        Returns:
            Dictionary mapping M-CSA ID to entry metadata

        Raises:
            FileNotFoundError: If no files match the pattern
            json.JSONDecodeError: If JSON is invalid
        """
        entry_info: Dict[int, dict] = {}
        for chunk in natsorted(glob(jsons_wildcard)):
            try:
                with open(chunk, 'r') as f:
                    for result in json.load(f)['results']:
                        entry_info[result['mcsa_id']] = result
            except (IOError, json.JSONDecodeError) as e:
                logger.warning(f"Error loading {chunk}: {e}")
        return entry_info

    @staticmethod
    def _get_residue_info(json_file: str) -> Dict[int, list]:
        """Load catalytic residue data from M-CSA JSON file.

        Parses M-CSA catalytic residue information including both reference
        and homologous residues, organized by entry ID.

        Args:
            json_file: Path to catalytic_residues_homologues.json

        Returns:
            Dictionary mapping M-CSA ID to list of residue records

        Raises:
            FileNotFoundError: If file not found
            json.JSONDecodeError: If JSON is invalid
        """
        json_residues: Dict[int, list] = defaultdict(list)
        try:
            with open(json_file, 'r') as f:
                for json_res in json.load(f):
                    json_residues[json_res['mcsa_id']].append(json_res)
        except (KeyError, IOError) as e:
            logger.error(f"Error loading residue info from {json_file}: {e}")
            return {}
        return json_residues

    @staticmethod
    def _get_pdb_assemblies(residue_info: Dict[int, list]) -> Dict[str, int]:
        """Map each PDB ID to the assembly number M-CSA annotates it in."""
        assemblies: Dict[str, int] = {}
        for json_residues in residue_info.values():
            for res in json_residues:
                for res_chain in res['residue_chains']:
                    assemblies.setdefault(res_chain['pdb_id'], res_chain.get('assembly') or 1)
        return assemblies

    def _get_cif_path(self, pdb_id: str) -> Optional[str]:
        """Path to the assembly mmCIF for a PDB entry, downloading it if missing.

        Returns None (with a warning) if the file is absent and cannot be fetched.
        """
        assemblies_dir = datasets.ASSEMBLIES_DIR
        matches = glob(os.path.join(assemblies_dir, f'{pdb_id}-assembly-*.cif'))
        if matches:
            return matches[0]
        path = fetch_assembly(pdb_id, self.pdb_assemblies.get(pdb_id, 1), assemblies_dir)
        if path is None:
            logger.warning(f'Could not find or download mmCIF file for {pdb_id}')
            return None
        return str(path)


__all__ = ['Mcsa']

