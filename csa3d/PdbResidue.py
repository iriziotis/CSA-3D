"""PDB catalytic residue representation for M-CSA."""

import logging
import os
import warnings
from typing import Optional, List, Dict, Any, Tuple, Generator, Union
import numpy as np
from copy import copy
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from Bio.PDB.Polypeptide import is_aa
from rdkit.Chem.rdmolfiles import MolFromMolFile
from rdkit.Chem import AllChem
from rdkit import DataStructs

from . import datasets
from .residue_definitions import (
    AA_3TO1,
    STANDARD_RESIDUES,
    NUCLEIC,
    EQUIVALENT_RESIDUES,
    EQUIVALENT_ATOMS,
    RESIDUE_DEFINITIONS,
)


logger = logging.getLogger(__name__)


class PdbResidue:
    """Catalytic residue in a PDB structure.

    Represents a single catalytic residue from the M-CSA database with 3D coordinates,
    structural information from PDB, and functional annotations. Each residue contains
    a reference to the parent site and optional Biopython Residue object with atomic details.

    Attributes:
        mcsa_id: Mechanism and Catalytic Site Atlas entry identifier
        pdb_id: PDB structure identifier (4-letter code, e.g., "1ABC")
        resname: Residue name in 3-letter code (e.g., "ALA", "SER"). Empty string for gaps.
        resid: Internal sequential residue number used in file formats
        auth_resid: Author-assigned residue number (e.g., displayed in PDB viewers)
        chain: Chain identifier (e.g., "A", "AB" for assembly chains)
        alt_chain: Alternative chain ID from the source file
        funclocs: List of functional location descriptions (e.g., ["Side chain"])
        is_reference: True if this residue is from the M-CSA reference structure
        is_standard: True if one of the 20 standard amino acids
        chiral_id: Identifier for stereoisomer when multiple catalytic roles exist
        structure: Biopython Residue object with full atomic coordinates and details
        reference_residue: Parent PdbResidue if this is a homologue
        parent_site: Parent PdbSite object
    """
    """Catalytic residue in a PDB structure.

    Represents a single catalytic residue with 3D coordinates, structural
    information, and functional annotations from M-CSA.

    Attributes:
        parent_site: Parent PdbSite object
        mcsa_id: M-CSA entry ID
        pdb_id: PDB identifier
        resname: Residue name (3-letter code)
        resid: Internal residue number
        auth_resid: Author-provided residue number
        chain: Chain identifier
        structure: Biopython Residue object
        reference_residue: Reference residue for comparison
    """

    def __init__(
        self,
        mcsa_id: Optional[int] = None,
        pdb_id: Optional[str] = None,
        resname: str = '',
        resid: Optional[int] = None,
        auth_resid: Optional[int] = None,
        chain: str = '',
        alt_chain: str = '',
        ptm_resname: str = '',
        funclocs: Optional[List[str]] = None,
        domain_name: str = '',
        domain_cath_id: str = '',
        is_reference: bool = False,
        chiral_id: Optional[int] = None,
        dummy_structure: bool = False,
    ) -> None:
        """Initialize a PdbResidue object.

        Args:
            mcsa_id: M-CSA entry identifier
            pdb_id: PDB structure ID
            resname: Residue 3-letter code
            resid: Internal residue number
            auth_resid: Author residue number
            chain: Chain identifier
            alt_chain: Alternative chain ID
            ptm_resname: Post-translational modification name if applicable
            funclocs: Functional location descriptions
            domain_name: Structural domain name
            domain_cath_id: CATH structural classification
            is_reference: Flag for M-CSA reference structure
            chiral_id: Stereoisomer identifier
            dummy_structure: If True, create a dummy atom for visualization
        """
        self.parent_site: Optional[Any] = None
        self.mcsa_id: Optional[int] = mcsa_id
        self.pdb_id: Optional[str] = pdb_id
        self.resname: str = resname
        self.resid: Optional[int] = resid
        self.auth_resid: Optional[int] = auth_resid
        self.chain: str = chain
        self.alt_chain: str = alt_chain
        self.ptm_resname: str = ptm_resname
        self.funclocs: List[str] = funclocs if funclocs else []
        self.domain_name: str = domain_name
        self.domain_cath_id: str = domain_cath_id
        self.is_reference: bool = is_reference
        self.chiral_id: Optional[int] = chiral_id
        self.reference_residue: Optional['PdbResidue'] = None
        self.structure: Optional[Residue] = None
        self.dummy_structure: Optional[Residue] = None
        if dummy_structure:
            self.add_dummy_structure()

    def __str__(self) -> str:
        """Return one-letter code representation."""
        return AA_3TO1.get(self.resname, 'X')

    def __eq__(self, other: Any) -> bool:
        """Check if residues share same information."""
        if not isinstance(other, PdbResidue):
            return False
        return (
            self.resname == other.resname
            and self.auth_resid == other.auth_resid
            and self.chain == other.chain
            and self.chiral_id == other.chiral_id
        )

    def __sub__(self, other: 'PdbResidue') -> float:
        """Return approximate distance to other residue."""
        return self.get_distance(other)

    # Private methods

    def _avg_distance(self, other: 'PdbResidue') -> float:
        """Average distance between two residues."""
        dists = []
        for i in self.structure.get_atoms():
            for j in other.structure.get_atoms():
                dists.append(i - j)
        return np.mean(np.array(dists))

    def _min_distance(self, other: 'PdbResidue') -> float:
        """Minimum distance between two residues."""
        min_dist = 999.0
        for i in self.structure.get_atoms():
            for j in other.structure.get_atoms():
                dist = i - j
                if dist < min_dist:
                    min_dist = dist
        return min_dist

    def _com_distance(
        self, other: 'PdbResidue', geometric: bool = True
    ) -> float:
        """Center of mass or geometry distance.

        Args:
            other: Residue to compare
            geometric: Use center of geometry (True) vs center of mass (False)

        Returns:
            Distance in Angstroms
        """
        if self.structure is None or other.structure is None:
            return np.nan
        self_com = self.structure.center_of_mass(geometric)
        other_com = other.structure.center_of_mass(geometric)
        return float(np.linalg.norm(self_com - other_com))

    # Public methods

    def copy(self, include_structure: bool = False) -> 'PdbResidue':
        """Create a shallow or deep copy of this residue.
        
        Args:
            include_structure: If True, also copy the Biopython Residue object
                with all atomic coordinates. If False (default), structure is set
                to None. Copying the structure is slower but needed for independent
                coordinate manipulation.
        
        Returns:
            PdbResidue: A copy of this residue with the same attributes.
        """
        res = copy(self)
        if include_structure:
            try:
                res.structure = self.structure.copy()
                res.structure.set_parent(self.structure.get_parent())
            except AttributeError:
                res.structure = None
        return res

    def add_structure(self, structure: 'Union[Residue, Structure]') -> bool:
        """Map this residue to atomic coordinates from a PDB structure.
        
        Finds and stores the corresponding Biopython Residue object that contains
        the atomic coordinates and connectivity information. If given a Structure
        object, automatically finds the residue with matching chain and position.
        
        Args:
            structure: Either a Biopython Residue object that matches this residue,
                or a Biopython Structure object from which the matching residue
                will be extracted using chain ID and author residue number.
        
        Returns:
            True if structure was successfully added or if this is a gap residue
            (no alignment). False if the residue cannot be found in the given
            structure or if coordinates cannot be extracted.
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

    def add_dummy_structure(self) -> bool:
        """Adds a dummy atom of zero coordinates to mark a gap in visualisation
        software"""
        dummy_atom = Atom('DUM', np.zeros(3), 0, 1, ' ', 'DUM', -999)
        dummy_residue = Residue((' ', -1 * self.chiral_id, ' '), 'DUM', '?')
        dummy_residue.add(dummy_atom)
        dummy_chain = Chain('?')
        dummy_chain.add(dummy_residue)
        self.dummy_structure = dummy_residue
        return True

    def get_distance(self, other: 'PdbResidue', kind: str = 'com') -> float:
        """Calculate distance to another residue.
        
        Args:
            other: The residue to compare distance to.
            kind: Type of distance metric:
                - 'com': Center of mass (default)
                - 'cog': Center of geometry
                - 'min': Minimum atom-to-atom distance
                - 'avg': Average atom-to-atom distance
                - 'ca': C-alpha distance (falls back to 'com' if missing)
        
        Returns:
            float: Distance in Ångströms, or np.nan if structures are unavailable.
        """
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

    def is_equivalent(self, other: 'PdbResidue', by_chiral_id: bool = True, by_chain: bool = False) -> bool:
        """Check if residues share the same pdb_id, chiral_id, name, resid
        and auth_resid"""
        basic = self.reference_residue == other.reference_residue
        chiral_ids = self.chiral_id == other.chiral_id
        chains = self.chain == other.chain

        if by_chiral_id:
            if by_chain:
                return basic and chiral_ids and chains
            return basic and chiral_ids
        if by_chain:
            return basic and chains
        return basic

    def get_func_atoms(self, allow_symmetrics: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """Extract functional atoms for superposition and analysis.
        
        Retrieves atoms that are relevant for catalytic function (excluding
        hydrogens and non-functional side chain atoms). Used for superposition
        and atom correspondence calculations.
        
        Args:
            allow_symmetrics: If True (default), treat symmetric atoms and
                residue equivalents as identical for matching. Uses definitions
                from residue_definitions.EQUIVALENT_ATOMS. If False, distinguishes
                all atoms individually.
        
        Returns:
            tuple: (atoms, coords) where:
                - atoms: NumPy array of atom ID strings, format "RESNAME.ATOMNAME"
                - coords: NumPy array of atomic coordinates (N×3), in Ångströms
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
    def id(self) -> Tuple[Optional[str], str, str, Optional[int], Optional[int]]:
        """ID of a residue as a tuple, not including chiral_id. Might be ambiguous"""
        return self.pdb_id, self.resname, self.chain, self.resid, self.auth_resid

    @property
    def full_id(self) -> Tuple[Optional[str], str, str, Optional[int], Optional[int], Optional[int]]:
        """ID of a residue as a tuple, including chiral_id. Always unique"""
        return self.pdb_id, self.resname, self.chain, self.resid, self.auth_resid, self.chiral_id

    @property
    def corrected_auth_resid(self) -> Optional[int]:
        """Correction of auth_resid for PDB files where identical assembly chains are
        identified by incrementing the auth_resid by 1000,2000 etc."""
        if not self.auth_resid:
            return
        if self.auth_resid < 1000 or len(self.chain) > 1:
            return self.auth_resid
        return int(str(self.auth_resid)[1:]) + (ord(self.chain) - 65) * 1000

    @property
    def is_gap(self) -> bool:
        """Check if residue is empty (no alignment with reference)"""
        return self.resname == '' and self.chain == '' and \
               self.resid is None and self.auth_resid is None

    @property
    def is_conserved(self) -> bool:
        """Check if residue are conserved by comparing to the reference"""
        return self.resname == self.reference_residue.resname or self.is_reference

    @property
    def is_conservative_mutation(self) -> bool:
        """Checks if residue is functionally equivalent to its reference"""
        if self.resname in EQUIVALENT_RESIDUES:
            if self.reference_residue.resname in EQUIVALENT_RESIDUES[self.resname]:
                return True
        return False

    @property
    def is_standard(self) -> bool:
        """Checks if residue is one of the 20 standard ones"""
        return self.resname.upper() in STANDARD_RESIDUES and \
               'ptm' not in self.funclocs and \
               (self.reference_residue is None or self.reference_residue.ptm_resname == '')

    @property
    def has_double_funcloc(self) -> bool:
        """Checks if residue has two function locations"""
        return len(self.funclocs) > 1

    @property
    def has_main_chain_function(self) -> bool:
        """Checks if residue has main chain function"""
        for funcloc in self.funclocs:
            if 'main' in funcloc.lower():
                return True
        return False

    @classmethod
    def from_json(cls, residue: Dict[str, Any], chiral_id: Optional[int] = None) -> Generator['PdbResidue', None, None]:
        """Generate PdbResidue objects from M-CSA JSON data.
        
        Parses a catalytic residue entry from the M-CSA database and yields
        PdbResidue objects for each occurrence in different PDB structures.
        
        Args:
            residue: A single catalytic residue entry from M-CSA JSON file
                (from catalytic_residues_homologues.json)
            chiral_id: Stereoisomer identifier (0 or 1 for distinct stereo roles).
                Defaults to None.
        
        Yields:
            PdbResidue: One residue object for each PDB structure containing
                this catalytic residue position.
        """
        for pdb_res in residue['residue_chains']:
            mcsa_id = residue['mcsa_id']
            pdb_id = pdb_res['pdb_id']
            resname = pdb_res['code']
            resid = pdb_res['resid']
            auth_resid = pdb_res['auth_resid']
            is_reference = pdb_res['is_reference']
            chain = PdbResidue.transform_chain(pdb_res['assembly_chain_name']) if is_reference \
                else PdbResidue.transform_chain(pdb_res['chain_name'])
            alt_chain = pdb_res['chain_name'] if is_reference \
                else pdb_res['assembly_chain_name']
            funclocs = [residue['function_location_abv']]
            domain_name = pdb_res['domain_name']
            domain_cath_id = pdb_res['domain_cath_id']
            #TODO wait for fixed json
            ptm_resname = residue.get('ptm', '')

            yield cls(mcsa_id, pdb_id, resname, resid, auth_resid,
                      chain, alt_chain, ptm_resname, funclocs, domain_name,
                      domain_cath_id, is_reference, chiral_id)
        return

    @staticmethod
    def transform_chain(chain: str, to_dash: bool = False) -> str:
        """Transform chain identifier between PDB and assembly formats.
        
        Converts chain IDs between two formats:
        - Standard: Single or double letter (e.g., "A", "AB")
        - Dash format: Letter-number with dash (e.g., "A-2")
        
        Useful for converting between conventional PDB chains and biological
        assembly notation where chains are numbered.
        
        Args:
            chain: Chain ID to transform
            to_dash: If True, convert to dash format "X-n". If False (default),
                convert from dash format to letter format "XY".
        
        Returns:
            str: Transformed chain ID, or original if already in target format.
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
    """Ligand, cofactor, ion, or substrate near a catalytic site.
    
    Represents non-protein components (heteroatoms) found within 5 Ångströms of
    catalytic residues. Extends PdbResidue with chemistry-specific classification
    and scoring. Automatically classifies components as artifacts, ions, cofactors,
    or substrates based on similarity to reaction components and known annotations.
    
    Attributes:
        similarity: Chemical fingerprint similarity (0-1) to the enzyme's known
            reaction components. Higher values indicate likely substrates or
            cofactors. None if structure files unavailable.
        best_match: Name of the best-matching reaction component (if similarity
            was computed).
        centrality: Mean minimum distance to catalytic residues, normalized by
            average intra-residue distance. Indicates how central the component
            is to the active site geometry.
        type: Classification as one of: "Ion", "Artefact", "Substrate (polymer)",
            "Substrate (non-polymer)", "Co-factor (ion)", "Co-factor (non-ion)".
        is_metal: True if single-atom metal (Fe, Zn, Mg, etc.)
        is_ion: True if ion or metal
        is_cofactor: True if annotated as cofactor in PDBe or M-CSA
        is_artefact: True if likely crystallographic artefact (e.g., PEG, glycerol)
    """

    def __init__(
        self,
        mcsa_id: Optional[int] = None,
        pdb_id: Optional[str] = None,
        resname: str = '',
        resid: Optional[int] = None,
        chain: str = '',
        structure: Optional[Residue] = None,
        parent_site: Optional[Any] = None,
        calculate_scores: bool = True
    ) -> None:
        """Initialize a Het (ligand/cofactor) object.
        
        Args:
            mcsa_id: M-CSA entry ID
            pdb_id: PDB structure identifier
            resname: Heteroatom residue name (e.g., "ZN", "NAG", "ADP")
            resid: Residue number in structure
            chain: Chain identifier
            structure: Biopython Residue object with atomic coordinates
            parent_site: Parent PdbSite object for calculating centrality
            calculate_scores: If True (default), compute similarity and centrality
                scores using the parent site. Set False for faster creation if
                these values are not needed.
        """
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

    def __repr__(self) -> str:
        return self.resname

    # Alternative constructors 

    @classmethod
    def polymer(cls, reslist: List[Residue], mcsa_id: Optional[int] = None, pdb_id: Optional[str] = None, chain: str = '', parent_site: Optional[Any] = None) -> 'Het':
        """Create a Het object representing a protein or nucleic acid chain.
        
        Used for substrates that are macromolecules (whole protein chain or
        nucleic acid) rather than small molecules.
        
        Args:
            reslist: List of Biopython Residue objects comprising the polymer
            mcsa_id: M-CSA entry ID
            pdb_id: PDB structure identifier
            chain: Chain identifier
            parent_site: Parent PdbSite for scoring
        
        Returns:
            Het: Object representing the polymer with resname='*P*'
        """
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
    def is_polymer(self) -> Optional[bool]:
        """Check if component is a peptide or nucleic macromolecule"""
        if self.structure:
            return type(self.structure) == Chain
        return None

    @property
    def is_peptide(self) -> Optional[bool]:
        """Check if component comes from a polypeptide"""
        if self.structure:
            return self.is_polymer and all([is_aa(res.get_resname()) for res in self.structure.get_residues()])
        return None

    @property
    def is_nucleic(self) -> Optional[bool]:
        """Check if component is DNA or RNA"""
        return self.is_polymer and not self.is_peptide

    @property
    def is_metal(self) -> bool:
        """Check if component is single-atom metal"""
        return self.resname in datasets.METALS

    @property
    def is_reactive_nonmetal(self) -> bool:
        """Check if component is a metal or a reactive non metal"""
        return self.resname in datasets.REACTIVE_NONMETALS

    @property
    def is_noble_gas(self) -> bool:
        """Check if component is a noble gas atom"""
        return self.resname in datasets.NOBLE_GASES

    @property
    def is_ion(self) -> bool:
        """Check if component is an ion"""
        return self.is_metal or self.is_reactive_nonmetal or self.is_noble_gas

    @property
    def is_cofactor(self) -> bool:
        """Check if component is annotated as cofactor for this pdb in PDBe or M-CSA"""
        return (self.resname in datasets.PDBID_COFACTORS[self.pdb_id]
                or self.resname in datasets.MCSA_EC_COFACTORS[self.parent_site.ec])

    @property
    def is_artefact(self) -> bool:
        """Check if component is likely to be a crystallographic artefact"""
        if self.cannot_be_artefact == True:
            return False
        return (self.resname in datasets.CRYSTALLIZATION_HETS and (self.similarity is None or self.similarity < 0.3))

    @property
    def type(self) -> str:
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

    def get_reaction(self) -> List[str]:
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

    def get_centrality(self) -> Optional[float]:
        """Calculate how central this ligand is to the catalytic site geometry.
        
        Computes the mean minimum distance from this component to all catalytic
        residues, normalized by the typical intra-residue distance. This gives
        a unitless measure where lower values indicate components more tightly
        bound into the active site geometry.
        
        Returns:
            float: Centrality score (unitless), or None if parent_site is unavailable.
        """
        dists = []
        if not self.parent_site:
            return None
        for res in self.parent_site:
            dists.append(self.get_distance(res, kind='com'))
        mean_site_dist = np.nanmean(self.parent_site.get_distances(kind='com'))
        centrality = np.nanmean(np.array(dists)) / mean_site_dist
        return np.round(centrality, 2)

    def get_similarity(self) -> Tuple[Optional[float], Optional[str]]:
        """Calculates the best PARITY match with the cognate reaction components"""
        #TODO: Check the dataset first and then calculate parity.
        if self.structure and self.structure.get_id()[0] != ' ':
            try:
                return self.generate_parity()
            except Exception as e:
                return None, None
        return None, None

    def generate_parity(self) -> Tuple[Optional[float], Optional[str]]:
        """Calculate Tanimoto similarity score for cognate reaction components.

        Computes molecular fingerprint similarity (Tanimoto) between the bound
        molecule and each cognate reaction component, returning the best match
        score and its ID.
        """
        if self.is_polymer:
            return None, None
        bound = f'{datasets.HET_MOLS_DIR}/{self.resname}.sdf'
        max_score = 0.0
        match = None
        components = self.get_reaction()
        if not components:
            return None, None

        try:
            bound_mol = MolFromMolFile(bound)
            if bound_mol is None or bound_mol.GetNumAtoms() == 0:
                return None, None
            bound_fp = AllChem.GetMorganFingerprintAsBitVect(bound_mol, 2, nBits=1024)
        except Exception:
            return None, None

        for cognate in components:
            try:
                cognate_mol = MolFromMolFile(cognate)
                if cognate_mol is None or cognate_mol.GetNumAtoms() == 0:
                    continue
                cognate_fp = AllChem.GetMorganFingerprintAsBitVect(cognate_mol, 2, nBits=1024)
                score = DataStructs.TanimotoSimilarity(bound_fp, cognate_fp)
            except Exception:
                continue

            if score >= max_score:
                max_score = score
                match = cognate.split('/')[-1].split('.')[0]

            # Check if cognate and bound component are single-atom molecules
            if self.is_polymer or (bound_mol.GetNumAtoms() == 1 and cognate_mol.GetNumAtoms() == 1):
                self.cannot_be_artefact = True

        return np.round(max_score, 2), match

    # Private methods

    def _get_mcsa_reaction(self) -> List[str]:
        """Gets component filenames of the parent M-CSA reaction."""
        components = []
        try:
            for component in self.parent_site.parent_entry.info['reaction']['compounds']:
                components.append('{}/ChEBI_{}.sdf'.format(datasets.CHEBI_MOLS_DIR, component['chebi_id']))
        except KeyError:
            pass
        return components

    def _get_kegg_reaction(self) -> List[str]:
        """Gets component filenames of the corresponding KEGG reaction."""
        return ['{}/KEGG_{}.sdf'.format(datasets.KEGG_MOLS_DIR, component)
                for component in datasets.KEGG_EC_REACTION.get(self.parent_site.ec, [])]

    def _get_rhea_reaction(self) -> List[str]:
        """Gets component filenames of the corresponding RHEA reaction."""
        return ['{}/ChEBI_{}.sdf'.format(datasets.CHEBI_MOLS_DIR, component)
                for component in datasets.RHEA_EC_REACTION.get(self.parent_site.ec, [])]


__all__ = ['PdbResidue', 'Het']
