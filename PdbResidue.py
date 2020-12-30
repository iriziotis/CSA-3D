# !/usr/bin/env python3

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Residue import Residue
from .residue_definitions import AA_3TO1, EQUIVALENT_RESIDUES
from ..templates.common import transform_chain


class PdbResidue:
    """M-CSA PDB residue"""

    def __init__(self, mcsa_id=None, pdb_id=None, resname='', resid=None,
                 auth_resid=None, chain='', funcloc=None, is_reference=False):
        self.mcsa_id = mcsa_id
        self.pdb_id = pdb_id
        self.resname = resname
        self.resid = resid
        self.auth_resid = auth_resid
        self.chain = chain
        self.funcloc = funcloc
        self.is_reference = is_reference
        self.reference_residue = None
        self.structure = None

    def __str__(self):
        return AA_3TO1[self.resname]

    def __eq__(self, other):
        return (self.resname == other.resname and
                self.auth_resid == other.auth_resid and
                self.chain == other.chain)

    def __sub__(self, other):
        return self.get_distance(other)

    def add_structure(self, path):
        """Map residue to Biopython residue object"""
        if path is None:
            return False
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure(self.pdb_id, path)
        try:
            structure = struct[0][self.chain][self.auth_resid]
            if structure.get_resname().capitalize() == self.resname:
                self.structure = structure
                return True
        except KeyError:
            return False

    def get_distance(self, other):
        """Get distance of two residues by measuring distance of two
        arbitrary atoms. Very approximate, use with care"""
        x = None
        y = None
        if self.structure is None or other.structure is None:
            return
        for atom in self.structure.get_atoms():
            x = atom
            break
        for atom in other.structure.get_atoms():
            y = atom
            break
        return x - y

    def is_equivalent(self, other):
        """Check if residues share the same pdb_id, name, resid
        and auth_resid"""
        return (self.pdb_id == other.pdb_id and
                self.resname == other.resname and
                self.auth_resid == other.auth_resid and
                self.resid == other.resid)

    def is_gap(self):
        """Check if residue is empty (no alignment with reference)"""
        return self.structure is None and self.mcsa_id is None and \
            self.pdb_id is None and self.resid is None

    def get_nearest_equivalent(self, other, reslist):
        """Finds the closest equivalent of another residue in a list"""
        min_dist = 999
        result = None
        for res in reslist:
            if res.structure is None:
                continue
            if res.is_equivalent(other):
                dist = self.get_distance(res)
                if dist < min_dist:
                    result = res
                    min_dist = dist
        return result

    @property
    def id(self):
        return self.pdb_id, self.resname, self.chain, self.resid, self.auth_resid

    @property
    def is_conserved(self):
        """Check if residue are conserved by comparing to the reference"""
        return self.resname == self.reference_residue.resname or self.is_reference

    def is_conservative_mutation(self):
        """Checks if residue is functionally equivalent to its reference"""
        if self.resname in EQUIVALENT_RESIDUES:
            if self.reference_residue.resname in EQUIVALENT_RESIDUES[self.resname]:
                return True
        return False

    @classmethod
    def from_json(cls, residue):
        """Constructs a list of PdbResidue objects using information directly from
        the M-CSA homologues json file. Input is a top-level residue entry in the json.
        """
        try:
            for pdb_res in residue['residue_chains']:

                mcsa_id = residue['mcsa_id']
                pdb_id = pdb_res['pdb_id']
                resname = pdb_res['code']
                resid = pdb_res['resid']
                auth_resid = pdb_res['auth_resid']
                is_reference = pdb_res['is_reference']
                chain = pdb_res['assembly_chain_name'] if is_reference \
                    else pdb_res['chain_name']
                chain = transform_chain(chain, to_dash=False)
                assembly_chain = transform_chain(pdb_res['assembly_chain_name'], to_dash=False)
                funcloc = residue['function_location_abv']

                if not is_reference and assembly_chain != chain:
                    yield cls(mcsa_id, pdb_id, resname, resid, auth_resid,
                              assembly_chain, funcloc, is_reference)
                yield cls(mcsa_id, pdb_id, resname, resid, auth_resid,
                          chain, funcloc, is_reference)
        except KeyError:
            return


class Het(PdbResidue):

    def __init__(self, mcsa_id=None, pdb_id=None, resname='', resid=None,
                 auth_resid=None, chain='', funcloc=None, is_reference=None,
                 parity_score=None, centrality=None):
        super().__init__(mcsa_id, pdb_id, resname, resid, chain)
        self.parity_score = parity_score
        self.centrality = centrality



