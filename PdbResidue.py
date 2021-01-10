from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from .residue_definitions import AA_3TO1, EQUIVALENT_RESIDUES


class PdbResidue:
    """M-CSA PDB residue. Basic information is collected from the M-CSA
    API catalytic residue info .json. Residue structure is set independently
    as a Biopython Residue object"""

    def __init__(self, mcsa_id=None, pdb_id=None, resname='', resid=None,
                 auth_resid=None, chain='', funcloc='', is_reference=False):
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
        """Returns one-letter representation"""
        return AA_3TO1[self.resname]

    def __eq__(self, other):
        """Checks if residues share the same information"""
        return (self.resname == other.resname and
                self.auth_resid == other.auth_resid and
                self.chain == other.chain)

    def __sub__(self, other):
        """Returns an approximate distance of self to other PdbResidue"""
        return self.get_distance(other)

    def copy(self):
        """Returns a copy of the object without the structure mapping"""
        res = PdbResidue(self.mcsa_id, self.pdb_id, self.resname, self.resid,
                         self.auth_resid, self.chain, self.funcloc, self.is_reference)
        res.reference_residue = self.reference_residue
        return res

    def add_structure(self, structure):
        """
        Map residue to Biopython residue object

        Args:
             structure: Biopython Structure or Residue object. If
                        object is Structure, it automatically gets the
                        Residue that has the same chain and auth_resid
        Returns:
            True if structure is added successfully or the PdbResidue is a gap
            (non-aligned, empty residue), False if residue cannot be found in the
            given structure.
        """
        if self.is_gap:
            return True
        if type(structure) == Residue:
            if structure.get_resname().capitalize() == self.resname and \
               structure.get_parent().get_id() == self.chain and \
               structure.get_id()[1] == self.auth_resid:
                self.structure = structure
                return True
        if type(structure) == Structure:
            try:
                residue = structure[0][self.chain][self.auth_resid]
                if residue.get_resname().capitalize() == self.resname:
                    self.structure = residue
                    return True
            except KeyError:
                return False
        else:
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
                #self.auth_resid == other.auth_resid and
                self.resid == other.resid)

    def get_nearest_equivalent(self, other, reslist):
        """
        Finds the closest equivalent of another residue in a list.
        Equivalent = same pdb_id, resid, auth_resid and resname,
                     different chain

        Args:
            other: The residue for which we want get the closest equivalent
                   in the list
            reslist: A list of PdbResidue objects
        Returns:
            The equivalent residue to other that is located closest to self
            in the 3D space.
        """
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
        """Unique ID of a residue as a tuple"""
        return self.pdb_id, self.resname, self.chain, self.resid, self.auth_resid

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
                funcloc = residue['function_location_abv']

                yield cls(mcsa_id, pdb_id, resname, resid, auth_resid,
                          chain, funcloc, is_reference)
        except KeyError:
            return

    @staticmethod
    def transform_chain(chain, to_dash=False):
        """
        Transforms the chain field from X-n to XY and reverse.
        Useful for chain IDs found in PDB assembly structures.

        Args:
            chain: The chain ID to be transformed
            to_dash: if True, transform to X-n format. If false,
                     transform to XY
        Returns:
            Transformed ID if input was of different format
        """
        letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        if to_dash:
            if len(chain) == 2 and all(c.isalpha for c in chain):
                return "{}-{}".format(chain[0], letters.find(chain[1])+2)
            else:
                return chain
        else:
            if '-' in chain:
                return "{}{}".format(chain.split("-")[0], letters[int(chain.split("-")[1])-2])
            else:
                return chain


class Het(PdbResidue):
    """PdbResidue with two extra attributes.
        parity_score: chemical similarity to the cognate ligand (PARITY score,
                      from a pre-compiled dataset)
        centrality: the mean minimum distance of the het to the catalytic
    """
    def __init__(self, mcsa_id=None, pdb_id=None, resname='', resid=None,
                 chain='', parity_score=None, centrality=None):
        super().__init__(mcsa_id, pdb_id, resname, resid, chain)
        self.parity_score = parity_score
        self.centrality = centrality
