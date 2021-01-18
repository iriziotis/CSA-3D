# !/usr/bin/env python3

from .residue_definitions import AA_3TO1


class UniResidue:
    """M-CSA UniProt residue. Basic information is collected from the M-CSA
    API catalytic residue info .json"""

    def __init__(self, mcsa_id=None, uniprot_id=None, resname='', resid=None,
                 funcloc=None, is_reference=False, index=None):
        self.mcsa_id = mcsa_id
        self.uniprot_id = uniprot_id
        self.resname = resname
        self.resid = resid
        self.funcloc = funcloc
        self.is_reference = is_reference
        self.index = index
        self.reference_residue = None

    def __str__(self):
        return AA_3TO1[self.resname]

    def __eq__(self, other):
        """Checks if residues share the same information"""
        return (self.uniprot_id == other.uniprot_id and
                self.resname == other.resname and
                self.index == other.index)

    def copy(self):
        """Returns a copy of the object without the structure mapping"""
        res = UniResidue(self.mcsa_id, self.uniprot_id, self.resname, self.resid,
                         self.funcloc, self.is_reference, self.index)
        res.reference_residue = self.reference_residue
        return res

    def is_equivalent(self, other, by_index=True):
        """Check if residues share the same pdb_id, resname, resid and optionally index"""
        basic = self.pdb_id == other.pdb_id and self.resname == other.resname and \
                self.resid == other.resid 
        indeces = self.index == other.index

        if by_index:
            return basic and indeces
        return basic

    @property
    def id(self):
        """Unique ID of a residue as a tuple"""
        return self.uniprot_id, self.resname, self.resid, self.index

    @property
    def is_gap(self):
        """Check if residue is empty (no alignment with reference)"""
        return self.mcsa_id is None and self.pdb_id is None and self.resid is None

    @property
    def is_conserved(self):
        """Check if residue are conserved by comparing to the reference"""
        return self.resname == self.reference_residue.resname or self.is_reference

    @property
    def is_conservative_mutation(self):
        """Checks if residue is functionally equivalent to its reference"""
        if self.resname in EQUIVALENT_RESIDUES:
            if self.reference_residue.resname in EQUIVALENT_RESIDUES[self.resname]:
                return True
        return False

    @classmethod
    def from_json(cls, residue, index=None):
        """Constructs a list of UniResidue objects using information directly from
        the M-CSA homologues .json file. Input is a top-level residue entry in the json.
        """
        try:
            for uniprot_res in residue['residue_sequences']:
                mcsa_id = residue['mcsa_id']
                uniprot_id = uniprot_res['uniprot_id']
                resname = uniprot_res['code']
                resid = uniprot_res['resid']
                is_reference = uniprot_res['is_reference']
                funcloc = residue['function_location_abv']

                yield UniResidue(mcsa_id, uniprot_id, resname, resid, 
                                 funcloc, is_reference, index)
        except KeyError:
            return
