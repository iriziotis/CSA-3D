import json
from collections import defaultdict
from glob import glob
from .config import INFO_JSON, ASSEMBLIES_DIR
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
        self.info_json = INFO_JSON
        self.entries = dict()
        self.json_residues = self._get_residue_info(self.info_json)
        self.pdb_residues = defaultdict(lambda: defaultdict(list))
        self.uni_residues = defaultdict(lambda: defaultdict(list))
        self.ref_pdb_residues = defaultdict(list)
        self.ref_uni_residues = defaultdict(list)

    def __iter__(self):
        yield from self.entries.values()

    def build(self, entries=None, annotate=True, redundancy_cutoff=None, verbose=False):
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
        if not entries:
            entries = self.json_residues.keys()
        if type(entries) != list:
            entries = [entries]
        for entry in entries:
            if str(entry) in self.json_residues.keys():
                self._build_pdb_residues(entry)
                self._build_uniprot_residues(entry)
                self._build_pdb_sites(entry, annotate, redundancy_cutoff, verbose)
                self._build_uniprot_sites(entry)
            else:
                return False
        return True

    def add(self, entry):
        """Adds entry to Mcsa object"""
        if type(entry) != Entry:
            return False
        self.entries[entry.mcsa_id] = entry
        return True

    # Properties

    @property
    def nofentries(self):
        return len(self.entries.keys())

    # Private methods

    def _build_pdb_residues(self, entry):
        """Builds PdbResidue objects from using raw info found in the .json
        file from M-CSA API (catalytic_residues_homologues.json). Extracts
        all equivalent residues from identical chains in biological assemblies,
        and creates individual instances for them"""
        reference_residue = None
        seen = set()
        for chiral_id, json_res in enumerate(self.json_residues[entry]):
            for residue in PdbResidue.from_json(json_res, chiral_id):
                if residue.is_reference:
                    reference_residue = residue
                    self.ref_pdb_residues[residue.mcsa_id].append(residue)
                else:
                    self.pdb_residues[residue.mcsa_id][residue.pdb_id].append(residue)
                residue.reference_residue = reference_residue

    def _build_uniprot_residues(self, entry):
        """Builds UniResidue objects from using raw info found in the .json
        file from M-CSA API (catalytic_residues_homologues.json)"""
        reference_residue = None
        for chiral_id, json_res in enumerate(self.json_residues[entry]):
            for residue in UniResidue.from_json(json_res, chiral_id):
                if residue.is_reference:
                    reference_residue = residue
                    self.ref_uni_residues[residue.mcsa_id].append(residue)
                else:
                    self.uni_residues[residue.mcsa_id][residue.uniprot_id].append(residue)
                residue.reference_residue = reference_residue

    def _build_pdb_sites(self, entry, annotate=True, redundancy_cutoff=None, verbose=False):
        """
        Builds PdbSite objects from PdbResidue lists using distance
        criteria and adds them to Entry objects.

        Args:
            annotate: To extract annotations (experimental method, resolution,
                      organism, host, EC etc.) for each catalytic at a slight
                      cost in execution time.
            verbose: Output the active site under process
        """
        self.entries[entry] = Entry(entry)
        try:
            ref_pdb_id = self.ref_pdb_residues[entry][0].pdb_id
        except (IndexError, KeyError):
            print('No reference PDB structure')
            return
        reference_site = PdbSite.build_reference(self.ref_pdb_residues[entry],
                                                 self._get_cif_path(ref_pdb_id), annotate)
        self.entries[entry].add(reference_site)
        for pdb_id, pdb in self.pdb_residues[entry].items():
            for site in PdbSite.build_all(pdb, reference_site, self._get_cif_path(pdb_id),
                                          annotate, redundancy_cutoff):
                if verbose:
                    print(site.id, site)
                self.entries[entry].add(site)
        # TODO check cases with two pdb references

    def _build_uniprot_sites(self, entry):
        """Builds UniSite objects from UniResidue lists and adds them
        to Entry objects"""
        if entry not in self.entries:
            self.entries[entry] = Entry(entry)
        reference_site = UniSite.build(self.ref_uni_residues[entry])
        self.entries[entry].add(reference_site)
        for uniprot_id, uniprot in self.uni_residues[entry].items():
            site = UniSite.build(uniprot, reference_site)
            self.entries[entry].add(site)
        # TODO see cases like mcsa 212 - Multiple uniprot reference seqs

    # Static methods

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
        return
