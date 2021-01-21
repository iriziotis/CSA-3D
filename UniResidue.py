# !/usr/bin/env python3

from .residue_definitions import AA_3TO1, STANDARD_RESIDUES, EQUIVALENT_RESIDUES


class UniResidue:
    """M-CSA UniProt residue. Basic information is collected from the M-CSA
    API catalytic residue info .json"""

    def __init__(self, mcsa_id=None, uniprot_id=None, resname='', resid=None,
                 funclocs=None, is_reference=False, chiral_id=None):
        self.mcsa_id = mcsa_id
        self.uniprot_id = uniprot_id
        self.resname = resname
        self.resid = resid
        self.funclocs = funclocs if funclocs else []
        self.is_reference = is_reference
        self.chiral_id = chiral_id
        self.reference_residue = None

    def __str__(self):
        return AA_3TO1[self.resname]

    def __eq__(self, other):
        """Checks if residues share the same information"""
        return (self.uniprot_id == other.uniprot_id and
                self.resname == other.resname and
                self.chiral_id == other.chiral_id)

    def copy(self):
        """Returns a copy of the object without the structure mapping"""
        res = UniResidue(self.mcsa_id, self.uniprot_id, self.resname, self.resid,
                         self.funclocs, self.is_reference, self.chiral_id)
        res.reference_residue = self.reference_residue
        return res

    def is_equivalent(self, other, by_chiral_id=True):
        """Check if residues share the same uniprot_id, resname, resid and optionally chiral_id"""
        basic = self.uniprot_id == other.uniprot_id and self.resname == other.resname and \
                self.resid == other.resid 
        chiral_ids = self.chiral_id == other.chiral_id

        if by_chiral_id:
            return basic and chiral_ids
        return basic

    @property
    def id(self):
        """ID of a residue as a tuple, not including chiral_id. Might be ambiguous"""
        return self.uniprot_id, self.resname, self.resid

    @property
    def full_id(self):
        """ID of a residue as a tuple, including chiral_id. Always unique"""
        return self.uniprot_id, self.resname, self.resid, self.chiral_id

    @property
    def is_gap(self):
        """Check if residue is empty (no alignment with reference)"""
        return self.resname == '' and self.resid is None

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

    @property
    def is_standard(self):
        """Checks if residue is one of the 20 standard ones"""
        return self.resname.upper() in STANDARD_RESIDUES

    @property
    def has_double_funcloc(self):
        """Checks if residue has two function locations"""
        return len(self.funclocs) > 1

    @property
    def has_main_chain_function(self):
        """Checks if residue has main chain function"""
        for funcloc in self.funclocs:
            if 'main' in funcloc.lower():
                return True
        return False

    @classmethod
    def from_json(cls, residue, chiral_id=None):
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
                funclocs = [residue['function_location_abv']]

                yield UniResidue(mcsa_id, uniprot_id, resname, resid, 
                                 funclocs, is_reference, chiral_id)
        except KeyError:
            return
