import warnings
import numpy as np
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from .residue_definitions import AA_3TO1, STANDARD_RESIDUES, EQUIVALENT_RESIDUES, EQUIVALENT_ATOMS, RESIDUE_DEFINITIONS
from .config import PDBID_COFACTORS, METAL_COFACTORS, CRYSTALLIZATION_HETS
from copy import copy


class PdbResidue:
    """M-CSA PDB residue. Basic information is collected from the M-CSA
    API catalytic residue info .json. Residue structure is set independently
    as a Biopython Residue object"""

    def __init__(self, mcsa_id=None, pdb_id=None, resname='', resid=None,
                 auth_resid=None, chain='', funclocs=None, is_reference=False, chiral_id=None,
                 dummy_structure=False):
        self.mcsa_id = mcsa_id
        self.pdb_id = pdb_id
        self.resname = resname
        self.resid = resid
        self.auth_resid = auth_resid
        self.chain = chain
        self.funclocs = funclocs if funclocs else []
        self.is_reference = is_reference
        self.chiral_id = chiral_id
        self.reference_residue = None
        self.structure = None
        self.dummy_structure = None
        if dummy_structure:
            self.add_dummy_structure()

    def __str__(self):
        """Returns one-letter representation"""
        return AA_3TO1[self.resname]

    def __eq__(self, other):
        """Checks if residues share the same information"""
        return (self.resname == other.resname and
                self.auth_resid == other.auth_resid and
                self.chain == other.chain and
                self.chiral_id == other.chiral_id)

    def __sub__(self, other):
        """Returns an approximate distance of self to other PdbResidue"""
        return self.get_distance(other)

    def copy(self, include_structure=False):
        """Returns a copy of the residue. If include_structure is False,
        then structure is not copied."""
        res = copy(self)
        if include_structure:
            try:
                res.structure = self.structure.copy()
                # This is to override an obscure biopython bug on disordered atoms
                # TODO Remove when they update Biopython
                res.structure.child_list = []
                res.structure.child_dict = {}
                for atom in self.structure.get_unpacked_list():
                    if atom.get_altloc() in (' ', 'A'):
                        res.structure.add(atom.copy())
                res.structure.set_parent(self.structure.get_parent())
            except AttributeError:
                res.structure = None
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
                    self.structure = residue.copy()
                    self.structure.set_parent(residue.get_parent())
                    return True
            except KeyError:
                if self.is_reference:
                    warnings.warn('Could not add residue {} structure'.format(self.id), RuntimeWarning)
                return False
        else:
            return False

    def add_dummy_structure(self):
        """Adds a dummy atom of zero coordinates to mark a gap in visualisation
        software"""
        dummy_atom = Atom('DUM', np.zeros(3), 0, 1, ' ', 'DUM', -999)
        dummy_residue = Residue((' ', -1 * self.chiral_id, ' '), 'DUM', '?')
        dummy_residue.add(dummy_atom)
        dummy_chain = Chain('?')
        dummy_chain.add(dummy_residue)
        self.dummy_structure = dummy_residue
        return True

    def get_distance(self, other, minimum=True):
        """Get distance of two residues. If minimum=True, minimum distance
        is calculated. If False, distance of CAs is returned"""
        # TODO make it use either CAs or minimum distance
        if self.structure is None or other.structure is None:
            return np.nan
        if minimum:
            min_dist = 999
            for i in self.structure.get_atoms():
                for j in other.structure.get_atoms():
                    dist = i - j
                    if dist < min_dist:
                        min_dist = dist
            return min_dist
        else:
            try:
                return self.structure['CA'] - other.structure['CA']
            except KeyError:
                return self.get_distance(other, minimum=True)

    def is_equivalent(self, other, by_chiral_id=True, by_chain=False):
        """Check if residues share the same pdb_id, chiral_id, name, resid
        and auth_resid"""
        basic = self.pdb_id == other.pdb_id and self.resname == other.resname and \
                (self.resid == other.resid or self.auth_resid == other.auth_resid)
        chiral_ids = self.chiral_id == other.chiral_id
        chains = self.chain == other.chain

        if by_chiral_id:
            if by_chain:
                return basic and chiral_ids and chains
            return basic and chiral_ids
        if by_chain:
            return basic and chains
        return basic

    def get_func_atoms(self, allow_symmetrics=True):
        """Gets atoms and coordinates for superposition and atom reordering
        calculations

        Args:
            allow_symmetrics: If True, equivalent residues and atoms
                              get the same id string, according to the
                              definitions in residue_definitions.py
                              (EQUIVALENT_ATOMS)
        Returns:
            atoms: A NumPy array of atom identifier strings of type
                   'N.RES.AT' where N is the residue serial number
                   in the .pdb file (consistent among all sites),
                   RES is the residue name and AT is the atom name
            coords: A NumPy array of the atomic coordinates
        """
        atoms = []
        coords = []
        if not self.structure:
            return np.array(atoms), np.array(coords)
        for atom in self.structure:
            resname = self.resname.upper()
            if allow_symmetrics:
                if self.has_main_chain_function:
                    resname = 'ANY'
                if not self.is_standard:
                    resname = 'PTM'
            atmid = '{}.{}'.format(resname, atom.name)
            if atmid in RESIDUE_DEFINITIONS:
                if allow_symmetrics:
                    if atmid in EQUIVALENT_ATOMS:
                        atmid = EQUIVALENT_ATOMS[atmid]
                atoms.append(atmid)
                coords.append(atom.get_coord())
        try:
            atoms = np.array(atoms, dtype=object)
            coords = np.stack(coords, axis=0)
        except ValueError:
            return None, None
        return atoms, coords

    @property
    def id(self):
        """ID of a residue as a tuple, not including chiral_id. Might be ambiguous"""
        return self.pdb_id, self.resname, self.chain, self.resid, self.auth_resid

    @property
    def full_id(self):
        """ID of a residue as a tuple, including chiral_id. Always unique"""
        return self.pdb_id, self.resname, self.chain, self.resid, self.auth_resid, self.chiral_id

    @property
    def corrected_auth_resid(self):
        """Correction of auth_resid for PDB files where identical assembly chains are
        identified by incrementing the auth_resid by 1000,2000 etc."""
        if not self.auth_resid:
            return
        if self.auth_resid < 1000 or len(self.chain) > 1:
            return self.auth_resid
        return int(str(self.auth_resid)[1:]) + (ord(self.chain) - 65) * 1000

    @property
    def is_gap(self):
        """Check if residue is empty (no alignment with reference)"""
        return self.resname == '' and self.chain == '' and \
               self.resid is None and self.auth_resid is None

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
                funclocs = [residue['function_location_abv']]

                yield cls(mcsa_id, pdb_id, resname, resid, auth_resid,
                          chain, funclocs, is_reference, chiral_id)
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
                return "{}-{}".format(chain[0], letters.find(chain[1]) + 2)
            else:
                return chain
        else:
            if '-' in chain:
                return "{}{}".format(chain.split("-")[0], letters[int(chain.split("-")[1]) - 2])
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
        super().__init__(mcsa_id=mcsa_id, pdb_id=pdb_id, resname=resname, resid=resid, chain=chain)
        self.parity_score = parity_score
        self.centrality = centrality

    @property
    def is_artefact(self):
        """Check if component is known to be used in protein crystallization media"""
        return self.resname in CRYSTALLIZATION_HETS

    @property
    def is_peptide(self):
        """Check if component is from a peptide moiety"""
        if self.structure:
            return self.structure.get_id()[0] == ' '
        return

    @property
    def is_cofactor(self):
        """Check if component is annotated as cofactor for this pdb in PDBe"""
        return self.resname in PDBID_COFACTORS[self.pdb_id]

    @property
    def is_metal(self):
        """Check if component is a metal or a metallic cofactor"""
        return self.resname in METAL_COFACTORS

    @property
    def flag(self):
        """An identifier to tell if component is an artefact, a peptide, a cofactor or a metallic
        compound (in priority order)."""
        flag = 'H'
        if self.is_artefact:
            flag = 'A'
        if self.is_metal:
            flag = 'M'
        if self.is_peptide:
            flag = 'P'
        if self.is_cofactor:
            flag = 'C'
        if self.is_cofactor and self.is_metal:
            flag = 'CM'
        return flag
