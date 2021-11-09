import os
import json
import warnings
import subprocess
from collections import defaultdict
from glob import glob
from natsort import natsorted
from .config import CAT_RES_INFO, MCSA_ENTRY_INFO, ASSEMBLIES_DIR, PDBID_RESOLUTION
from .Entry import Entry
from .PdbSite import PdbSite
from .UniSite import UniSite
from .PdbResidue import PdbResidue
from .UniResidue import UniResidue


class Mcsa:
    """Top-level database Class. Contains Entry objects and methods
    to add or build entries using .json formatted data from the M-CSA
    API. Pathnames of relevant datasets are defined in config.py module.
    Definitions of residues and atoms for fitting methods are defined in
    the residue_definitions.py module."""

    def __init__(self):
        self.entries = dict()
        self.mcsa_entry_info = self._get_mcsa_entry_info(MCSA_ENTRY_INFO)
        self.catalytic_residue_info = self._get_residue_info(CAT_RES_INFO)
        self.pdb_residues = defaultdict(lambda: defaultdict(list))
        self.uni_residues = defaultdict(lambda: defaultdict(list))
        self.ref_pdb_residues = defaultdict(list)
        self.ref_uni_residues = defaultdict(list)

    def __iter__(self):
        yield from self.entries.values()

    def build(self, entry_ids=None, annotate=True, redundancy_cutoff=None, resolution_cutoff=None,  verbose=False, no_sites=False):
        """
        Builds M-CSA entries (Entry objects) containing PDB and UniProt
        catalytic sites (PdbSite and UniSite objects respectively

        Args:
            entries: A list of integer M-CSA IDs to specify which entries will
                     be built. If None, all entries will be built.
            annotate: To extract annotations (experimental method, resolution,
                      organism, host, EC etc.) for each catalytic at a slight
                      cost in execution time.
            verbose: Output the active site under process
        """
        if not entry_ids:
            entry_ids = self.catalytic_residue_info.keys()
        if type(entry_ids) != list:
            entry_ids = [entry_ids]
        for entry_id in entry_ids:
            if entry_id in self.catalytic_residue_info.keys():
                self.entries[entry_id] = Entry(entry_id)
                self.entries[entry_id].info = self.mcsa_entry_info[entry_id]
                self._build_pdb_residues(entry_id, resolution_cutoff)
                self._build_uniprot_residues(entry_id)

                # TODO check cases with two pdb references
                if len(set([r.pdb_id for r in self.ref_pdb_residues[entry_id]]))>1:
                    print('Has multiple PDB references. Not yet implemented')
                    return False
                # Temporary, until I implement treating of multiple references
                # TODO see cases like mcsa 212 - Multiple uniprot reference seqs
                if len(set([r.uniprot_id for r in self.ref_uni_residues[entry_id]]))>1:
                    print('Has multiple UniProt references. Not yet implemented')
                    return False
                if no_sites:
                    continue
                self._build_pdb_sites(entry_id, annotate, redundancy_cutoff, verbose)
                self._build_uniprot_sites(entry_id)
            else:
                return False
        return True

    def add(self, entry):
        """Adds entry to Mcsa object"""
        if type(entry) != Entry:
            return False
        self.entries[entry.mcsa_id] = entry
        return True

    def update_pdbs(self):
        """Updates the raw assembly directory from PDBe"""
        pdbs = set()
        for entry, json_res in self.catalytic_residue_info.items():
            for res in json_res:
                for res_chain in res['residue_chains']:
                    pdb_id = res_chain['pdb_id']
                    assembly = res_chain['assembly'] if res_chain['assembly'] else 1
                    pdbs.add((pdb_id, assembly))
        for pdb in pdbs:
            pdb_id, assembly = pdb[0], pdb[1]
            filename = '{}-assembly-{}.cif.gz'.format(pdb_id, assembly)
            if os.path.isfile(os.path.join(ASSEMBLIES_DIR, filename)):
                continue
            print('Getting PDB entry {}, assembly {}'.format(pdb_id, assembly))
            scp_path = "/nfs/services/pdbe/release-data/www-static-content/entry/{}/{}/{}".format(pdb_id[1:3], pdb_id, filename)
            cmd = "cp {} {}".format(scp_path, ASSEMBLIES_DIR)
            subprocess.call(cmd.split(" "))
        return

    # Properties

    @property 
    def nofentries(self):
        return len(self.entries.keys())

    # Private methods

    def _build_pdb_residues(self, entry_id, resolution_cutoff=None):
        """Builds PdbResidue objects from using raw info found in the .json
        file from M-CSA API (catalytic_residues_homologues.json). Extracts
        all equivalent residues from identical chains in biological assemblies,
        and creates individual instances for them"""
        reference_residue = None
        for chiral_id, json_res in enumerate(self.catalytic_residue_info[entry_id]):
            for residue in PdbResidue.from_json(json_res, chiral_id):
                # Resolution cutoff
                if resolution_cutoff and not residue.is_reference:
                    resolution = float(PDBID_RESOLUTION.get(residue.pdb_id, -1.0))
                    if resolution == -1.0 or resolution > resolution_cutoff:
                        continue
                if residue.is_reference:
                    reference_residue = residue
                    self.ref_pdb_residues[residue.mcsa_id].append(residue)
                else:
                    self.pdb_residues[residue.mcsa_id][residue.pdb_id].append(residue)
                residue.reference_residue = reference_residue

    def _build_uniprot_residues(self, entry_id):
        """Builds UniResidue objects from using raw info found in the .json
        file from M-CSA API (catalytic_residues_homologues.json)"""
        reference_residue = None
        for chiral_id, json_res in enumerate(self.catalytic_residue_info[entry_id]):
            for residue in UniResidue.from_json(json_res, chiral_id):
                if residue.is_reference:
                    reference_residue = residue
                    self.ref_uni_residues[residue.mcsa_id].append(residue)
                else:
                    self.uni_residues[residue.mcsa_id][residue.uniprot_id].append(residue)
                residue.reference_residue = reference_residue

    def _build_pdb_sites(self, entry_id, annotate=True, redundancy_cutoff=None, verbose=False):
        """
        Builds PdbSite objects from PdbResidue lists using distance
        criteria and adds them to Entry objects.

        Args:
            annotate: To extract annotations (experimental method, resolution,
                      organism, host, EC etc.) for each catalytic at a slight
                      cost in execution time.
            verbose: Output the active site under process
        """
        try:
            ref_pdb_id = self.ref_pdb_residues[entry_id][0].pdb_id
        except (IndexError, KeyError):
            print('No reference PDB structure')
            return
        reference_site = PdbSite.build_reference(self.ref_pdb_residues[entry_id], self.entries[entry_id],
                                                 self._get_cif_path(ref_pdb_id), annotate)
        self.entries[entry_id].add(reference_site)
        for pdb_id, pdb in self.pdb_residues[entry_id].items():
            for site in PdbSite.build_all(pdb, reference_site, self.entries[entry_id],
                                          self._get_cif_path(pdb_id), annotate, redundancy_cutoff):
                if verbose:
                    print(site.id, site)
                self.entries[entry_id].add(site)

    def _build_uniprot_sites(self, entry_id):
        """Builds UniSite objects from UniResidue lists and adds them
        to Entry objects"""
        reference_site = UniSite.build(self.ref_uni_residues[entry_id])
        reference_site.parent_entry = self.entries[entry_id]
        self.entries[entry_id].add(reference_site)
        for uniprot_id, uniprot in self.uni_residues[entry_id].items():
            site = UniSite.build(uniprot, reference_site)
            site.parent_entry = self.entries[entry_id]
            self.entries[entry_id].add(site)

    # Static methods

    @staticmethod
    def _get_mcsa_entry_info(jsons_wildcard):
        """
        Get all basic info for an entry from the Json
        """
        # Make dict with entry info
        entry_info = {}
        for chunk in natsorted(glob(jsons_wildcard)):
            with open(chunk, 'r') as f:
                for result in json.load(f)['results']:
                    entry_info[result['mcsa_id']] = result
        return entry_info

    @staticmethod
    def _get_residue_info(json_file):
        """Parses the raw .json catalytic residue info file and returns a
        dict of homologous residue lists for each M-CSA ID in .json format"""
        json_residues = defaultdict(list)
        try:
            with open(json_file, 'r') as f:
                for json_res in json.load(f):
                    json_residues[json_res['mcsa_id']].append(json_res)
        except KeyError:
            return
        return json_residues

    @staticmethod
    def _get_cif_path(pdb_id):
        """Returns the parent structure mmCIF file path"""
        path = '{0}/{1}-assembly-*.cif'.format(ASSEMBLIES_DIR, pdb_id)
        if len(glob(path)) > 0:
            return glob(path)[0]
        warnings.warn('Could not find mmCIF file for {}'.format(pdb_id), RuntimeWarning)
        return

