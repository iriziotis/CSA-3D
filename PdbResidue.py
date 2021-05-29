import os
import warnings
import numpy as np
import parity_core
from copy import copy
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import is_aa
from rdkit.Chem.rdmolfiles import MolFromMolFile
from .residue_definitions import AA_3TO1, STANDARD_RESIDUES, NUCLEIC, EQUIVALENT_RESIDUES, EQUIVALENT_ATOMS, RESIDUE_DEFINITIONS
from .config import METALS, REACTIVE_NONMETALS, NOBLE_GASES, PDBID_COFACTORS, MCSA_EC_COFACTORS,\
                    CRYSTALLIZATION_HETS, PDB2EC, KEGG_EC_REACTION, RHEA_EC_REACTION,\
                    HET_MOLS_DIR, CHEBI_MOLS_DIR, KEGG_MOLS_DIR


class PdbResidue:
    """M-CSA PDB residue. Basic information is collected from the M-CSA
    API catalytic residue info .json. Residue structure is set independently
    as a Biopython Residue object"""

    def __init__(self, mcsa_id=None, pdb_id=None, resname='', resid=None,
                 auth_resid=None, chain='', alt_chain='', funclocs=None, 
                 is_reference=False, chiral_id=None, dummy_structure=False):
        self.parent_site = None
        self.mcsa_id = mcsa_id
        self.pdb_id = pdb_id
        self.resname = resname
        self.resid = resid
        self.auth_resid = auth_resid
        self.chain = chain
        self.alt_chain = alt_chain
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

    # Private methods

    def _avg_distance(self, other):
        """Average distance between two residues"""
        dists = []
        for i in self.structure.get_atoms():
            for j in other.structure.get_atoms():
                dists.append(i - j)
        return np.mean(np.array(dists))

    def _min_distance(self, other):
        """Minimum distance between two residues"""
        min_dist = 999
        for i in self.structure.get_atoms():
            for j in other.structure.get_atoms():
                dist = i - j
                if dist < min_dist:
                    min_dist = dist
        return min_dist

    def _com_distance(self, other, geometric=True):
        """Center of mass or geometry distance. If geometric=True, the
        distance is between the centers of geometry, otherwise it is between
        the centers of mass."""
        if self.structure is None or other.structure is None:
            return np.nan
        self_com = self.structure.center_of_mass(geometric)
        other_com = other.structure.center_of_mass(geometric)
        return np.linalg.norm(self_com-other_com)


    # Public methods

    def copy(self, include_structure=False):
        """Returns a copy of the residue. If include_structure is False,
        then structure is not copied."""
        res = copy(self)
        if include_structure:
            try:
                res.structure = self.structure.copy()
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
            if self.is_standard:
                try:
                    residue = structure[0][self.chain][self.auth_resid]
                except KeyError:
                    if self.is_reference:
                        warnings.warn('Could not add residue {} structure'.format(self.id), RuntimeWarning)
                    return False
            else:
                found = False
                for res in structure[0][self.chain].get_residues():
                    if res.get_id()[1] == self.auth_resid:
                        residue = res
                        found = True
                        break
                if not found and self.is_reference:
                    warnings.warn('Could not add residue {} structure'.format(self.id), RuntimeWarning)
                    return False
            if residue.get_resname().capitalize() == self.resname or not self.is_standard:
                self.structure = residue.copy()
                self.structure.set_parent(residue.get_parent())
                return True
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

    def get_distance(self, other, kind='com'):
        if self.structure is None or other.structure is None:
            return np.nan
        if kind == 'com':
            return self._com_distance(other, geometric=False)
        if kind == 'cog':
            return self._com_distance(other, geometric=True)
        if kind == 'min':
            return self._min_distance(other)
        if kind == 'avg':
            return self._avg_distance(other)
        if kind == 'ca':
            try:
                return self.structure['CA'] - other.structure['CA']
            except KeyError:
                return self.get_distance(other, kind='com')

    def is_equivalent(self, other, by_chiral_id=True, by_chain=False):
        """Check if residues share the same pdb_id, chiral_id, name, resid
        and auth_resid"""
        basic = self.pdb_id == other.pdb_id and \
                (self.resname == other.resname or not self.is_standard or not other.is_standard) and \
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
        return self.resname.upper() in STANDARD_RESIDUES and 'ptm' not in self.funclocs

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
                alt_chain = pdb_res['chain_name'] if is_reference \
                    else pdb_res['assembly_chain_name']
                funclocs = [residue['function_location_abv']]

                yield cls(mcsa_id, pdb_id, resname, resid, auth_resid,
                          chain, alt_chain, funclocs, is_reference, chiral_id)
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
                 chain='', structure=None, parent_site=None, calculate_scores=True):
        super().__init__(mcsa_id=mcsa_id, pdb_id=pdb_id, resname=resname, resid=resid, chain=chain)
        self.structure = None
        self.parent_site = parent_site
        self.cannot_be_artefact = False
        self.components_missing = False
        self.is_distal = False
        self.similarity, self.best_match = None, None
        self.centrality = None
        if structure:
            self.structure = structure.copy()
            self.structure.set_parent(structure.get_parent())
        if calculate_scores:
            self.similarity, self.best_match = self.get_similarity()
            self.centrality = self.get_centrality()

    def __repr__(self):
        return self.resname

    # Alternative constructors 

    @classmethod
    def polymer(cls, reslist, mcsa_id=None, pdb_id=None, chain='', parent_site=None):
        """Alternative constructor for polymers. Takes a residue list and returns
        a polymer ligand"""
        poly = cls(mcsa_id, pdb_id, resname='*P*', resid=None, chain=chain, 
                   structure=None, parent_site=parent_site, calculate_scores=False)
        poly.structure = Chain(chain)
        for res in reslist:
            if res.get_id() not in poly.structure:
                poly.structure.add(res.copy())
        poly.similarity, poly.best_match = poly.get_similarity()
        poly.centrality = poly.get_centrality()
        return poly

    # Properties

    @property
    def is_polymer(self):
        """Check if component is a peptide or nucleic macromolecule"""
        if self.structure:
            return type(self.structure) == Chain
        return

    @property
    def is_peptide(self):
        """Check if component comes from a polypeptide"""
        if self.structure:
            return self.is_polymer and all([is_aa(res.get_resname()) for res in self.structure.get_residues()])

    @property
    def is_nucleic(self):
        """Check if component is DNA or RNA"""
        return self.is_polymer and not self.is_peptide

    @property
    def is_metal(self):
        """Check if component is single-atom metal"""
        return self.resname in METALS

    @property
    def is_reactive_nonmetal(self):
        """Check if component is a metal or a reactive non metal"""
        return self.resname in REACTIVE_NONMETALS

    @property
    def is_noble_gas(self):
        """Check if component is a noble gas atom"""
        return self.resname in NOBLE_GASES

    @property
    def is_ion(self):
        """Check if component is an ion"""
        return self.is_metal or self.is_reactive_nonmetal or self.is_noble_gas

    @property
    def is_cofactor(self):
        """Check if component is annotated as cofactor for this pdb in PDBe or M-CSA"""
        return self.resname in PDBID_COFACTORS[self.pdb_id] or self.resname in MCSA_EC_COFACTORS[self.parent_site.ec]

    @property
    def is_artefact(self):
        """Check if component is likely to be a crystallographic artefact"""
        if self.cannot_be_artefact == True:
            return False
        return (self.resname in CRYSTALLIZATION_HETS and (self.similarity is None or self.similarity < 0.3))

    @property
    def type(self):
        """An identifier to tell if component is an artefact, a peptide, a cofactor or a metallic
        compound (in priority order)."""
        flag = 'Substrate (non-polymer)'
        if self.is_ion:
            flag = 'Ion'
        if self.is_artefact:
            flag = 'Artefact'
        if self.is_polymer:
            flag = 'Substrate (polymer)'
        if self.is_cofactor:
            if self.is_ion:
                flag = 'Co-factor (ion)'
            else:
                flag = 'Co-factor (non-ion)'
        return flag

    # Public methods

    def get_reaction(self):
        """Gets filenames of reactants and products either from M-CSA or KEGG, 
        depending on the EC match with the reference active site."""
        if self.parent_site.ec == self.parent_site.parent_entry.info['reaction']['ec']:
            components = self._get_mcsa_reaction()
        else:
            components = self._get_rhea_reaction()
            if not components:
                components = self._get_kegg_reaction()
        for component in reversed(components):
            if not os.path.exists(component):
                components.remove(component)
                self.components_missing = True
        return components

    def get_centrality(self):
        """Calculates an average distance to the catalytic residues of the site, 
        normalized by the average intra-residue distance of the site"""
        dists = [] 
        if not self.parent_site:
            return
        for res in self.parent_site:
            dists.append(self.get_distance(res, kind='com'))
        mean_site_dist = np.nanmean(self.parent_site.get_distances(kind='com'))
        centrality = np.nanmean(np.array(dists)) / mean_site_dist
        return np.round(centrality, 2)

    def get_similarity(self):
        """Calculates the best PARITY match with the cognate reaction components"""
        #TODO: Check the dataset first and then calculate parity.
        if self.structure and self.structure.get_id()[0] != ' ':
            try:
                return self.generate_parity()
            except Exception as e:
                return None, None
        return None, None
        
    def generate_parity(self):
        """Calculates parity score for each cognate reaction component and
        for the best match, it returns the score, the match ID, and the type
        (Reactant or Product)."""
        if self.is_polymer:
            return None, None
        bound = f'{HET_MOLS_DIR}/{self.resname}.sdf'
        max_score = 0.0
        match = None
        components = self.get_reaction()
        if not components:
            return None, None
        for cognate in components:
            try:
                score = parity_core.generate_parity(bound, cognate, quiet=True)
            except Exception as e:
                continue
            if score >= max_score:
                max_score = score
                match = cognate.split('/')[-1].split('.')[0]
            # Check if cognate and bound component are single-atom molecules
            if self.is_polymer or (MolFromMolFile(bound).GetNumAtoms() == 1 and MolFromMolFile(cognate).GetNumAtoms() == 1):
                self.cannot_be_artefact = True
        return np.round(max_score, 2), match

    # Private methods

    def _get_mcsa_reaction(self):
        """Gets component filenames of the parent M-CSA reaction."""
        components = []
        try:
            for component in self.parent_site.parent_entry.info['reaction']['compounds']:
                components.append('{}/ChEBI_{}.sdf'.format(CHEBI_MOLS_DIR, component['chebi_id']))
        except KeyError:
            pass
        return components

    def _get_kegg_reaction(self):
        """Gets component filenames of the corresponding KEGG reaction."""
        components = []
        try:
            for component in KEGG_EC_REACTION[self.parent_site.ec]:
                components.append('{}/KEGG_{}.mol'.format(KEGG_MOLS_DIR, component))
        except KeyError:
            pass
        return components

    def _get_rhea_reaction(self):
        """Gets component filenames of the corresponding KEGG reaction."""
        components = []
        try:
            for component in RHEA_EC_REACTION[self.parent_site.ec]:
                components.append('{}/ChEBI_{}.sdf'.format(CHEBI_MOLS_DIR, component))
        except KeyError:
            pass
        return components
