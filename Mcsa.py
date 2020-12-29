#!/usr/bin/env python3

import json
from collections import defaultdict
from glob import glob
from .config import INFO_JSON, WORKING_DIR, ASSEMBLIES_DIR
from .Entry import Entry
from .PdbSite import PdbSite
from .UniSite import UniSite
from .PdbResidue import PdbResidue
from .UniResidue import UniResidue


class Mcsa:

    def __init__(self):
        self.info_json = INFO_JSON
        self.nofentries = 0
        self.entries = dict()
        self.json_residues = self._get_residue_info(self.info_json)
        self.pdb_residues = defaultdict(lambda: defaultdict(list))
        self.uni_residues = defaultdict(lambda: defaultdict(list))
        self.ref_pdb_residues = defaultdict(list)
        self.ref_uni_residues = defaultdict(list)

    def __iter__(self):
        yield from self.entries.values()

    def build(self, entries=None):
        if not entries:
            entries = self.json_residues.keys()

        for entry in entries:
            self._build_pdb_residues(entry)
            self._build_uniprot_residues(entry)
            self._build_pdb_sites(entry)
            self._build_uniprot_sites(entry)
            
    def _build_pdb_residues(self, entry):
        reference_residue = None
        for json_res in self.json_residues[entry]:
            for residue in PdbResidue.from_json(json_res):
                if residue.is_reference:
                    reference_residue = residue
                    self.ref_pdb_residues[residue.mcsa_id].append(residue)
                else:
                    self.pdb_residues[residue.mcsa_id][residue.pdb_id].append(residue)
                residue.add_structure(self._get_catres_path(residue.mcsa_id, residue.pdb_id))
                residue.reference_residue = reference_residue

    def _build_uniprot_residues(self, entry):
        reference_residue = None
        for json_res in self.json_residues[entry]:
            for residue in UniResidue.from_json(json_res):
                if residue.is_reference:
                    reference_residue = residue
                    self.ref_uni_residues[residue.mcsa_id].append(residue)
                else:
                    self.uni_residues[residue.mcsa_id][residue.uniprot_id].append(residue)
                residue.reference_residue = reference_residue

    def _build_pdb_sites(self, entry):
        self.entries[entry] = Entry(entry)
        reference_site = PdbSite.build_reference(self.ref_pdb_residues[entry])
        self.entries[entry].add(reference_site)
        for pdb_id, pdb in self.pdb_residues[entry].items():
            for site in PdbSite.build_all(pdb, reference_site):
                cif = glob('{}/*{}*.cif'.format(ASSEMBLIES_DIR, site.pdb_id))[0]
                site.find_ligands(parent_structure=cif)
                self.entries[entry].add(site)
        # TODO check cases like mcsa 2, 1xxm, and cases with two pdb references

    def _build_uniprot_sites(self, entry):
        if entry not in self.entries:
            self.entries[entry] = Entry(entry)
        reference_site = UniSite.build(self.ref_uni_residues[entry])
        for uniprot_id, uniprot in self.uni_residues[entry].items():
            site = UniSite.build(uniprot, reference_site)
            self.entries[entry].add(site)
        # TODO see cases like mcsa 212 - Multiple uniprot reference seqs

    @staticmethod
    def _get_residue_info(json_file):
        json_residues = defaultdict(list)
        try:
            with open(json_file, 'r') as f:
                from .Entry import Entry
                for json_res in json.load(f):
                    json_residues[json_res['mcsa_id']].append(json_res)
        except KeyError:
            return
        return json_residues

    @staticmethod
    def _get_catres_path(mcsa_id, pdb_id):
        path = '{0}/results/entries/mcsa_{1}/mcsa_{1}.{2}.*.pdb'.format(
            WORKING_DIR, str(mcsa_id).zfill(4), pdb_id)
        if len(glob(path)) > 0:
            return glob(path)[0]
        return
