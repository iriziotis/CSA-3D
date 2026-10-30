# !/usr/bin/env python3
"""UniProt sequence residue representation for M-CSA catalytic sites."""

from typing import Optional, List, Tuple, Generator, Any, Dict

# Handle imports from both package and root-level contexts
try:
    from residue_definitions import AA_3TO1, STANDARD_RESIDUES, EQUIVALENT_RESIDUES
except ImportError:
    from .residue_definitions import AA_3TO1, STANDARD_RESIDUES, EQUIVALENT_RESIDUES


class UniResidue:
    """Catalytic residue in UniProt sequence.

    Represents a catalytic residue position in a UniProt protein sequence,
    with alignment to reference and functional annotations.

    Attributes:
        mcsa_id: M-CSA entry ID
        uniprot_id: UniProt accession
        resname: Residue name (3-letter code)
        resid: Position in sequence
        funclocs: Functional location descriptions
        is_reference: Whether this is the reference sequence
        chiral_id: Stereo identifier
        reference_residue: Reference residue for comparison
    """

    def __init__(
        self,
        mcsa_id: Optional[int] = None,
        uniprot_id: Optional[str] = None,
        resname: str = '',
        resid: Optional[int] = None,
        funclocs: Optional[List[str]] = None,
        is_reference: bool = False,
        chiral_id: Optional[int] = None,
    ) -> None:
        """Initialize UniResidue."""
        self.mcsa_id: Optional[int] = mcsa_id
        self.uniprot_id: Optional[str] = uniprot_id
        self.resname: str = resname
        self.resid: Optional[int] = resid
        self.funclocs: List[str] = funclocs if funclocs else []
        self.is_reference: bool = is_reference
        self.chiral_id: Optional[int] = chiral_id
        self.reference_residue: Optional['UniResidue'] = None
        self.mapped_pdbsites: List[Any] = []

    def __str__(self) -> str:
        """Return single-letter code."""
        return AA_3TO1.get(self.resname, 'X')

    def __eq__(self, other: Any) -> bool:
        """Check if residues share same information."""
        if not isinstance(other, UniResidue):
            return False
        return (
            self.uniprot_id == other.uniprot_id
            and self.resname == other.resname
            and self.chiral_id == other.chiral_id
        )

    def copy(self) -> 'UniResidue':
        """Return a copy without structure mapping."""
        res = UniResidue(
            self.mcsa_id,
            self.uniprot_id,
            self.resname,
            self.resid,
            self.funclocs.copy(),
            self.is_reference,
            self.chiral_id,
        )
        res.reference_residue = self.reference_residue
        return res

    def is_equivalent(self, other: 'UniResidue', by_chiral_id: bool = True) -> bool:
        """Check if residues are equivalent.

        Args:
            other: UniResidue to compare
            by_chiral_id: Include chiral ID in comparison. Defaults to True.

        Returns:
            True if equivalent
        """
        basic = (
            self.uniprot_id == other.uniprot_id
            and self.resname == other.resname
            and self.resid == other.resid
        )
        chiral_ids = self.chiral_id == other.chiral_id

        if by_chiral_id:
            return basic and chiral_ids
        return basic

    @property
    def id(self) -> Tuple[Optional[str], str, Optional[int]]:
        """Residue ID (ambiguous, may not include chiral_id)."""
        return self.uniprot_id, self.resname, self.resid

    @property
    def full_id(self) -> Tuple[Optional[str], str, Optional[int], Optional[int]]:
        """Unique residue ID including chiral_id."""
        return self.uniprot_id, self.resname, self.resid, self.chiral_id

    @property
    def is_gap(self) -> bool:
        """Check if residue is gap (no alignment)."""
        return self.resname == '' and self.resid is None

    @property
    def is_conserved(self) -> bool:
        """Check if residue is conserved vs reference."""
        if self.reference_residue is None:
            return self.is_reference
        return self.resname == self.reference_residue.resname or self.is_reference

    @property
    def is_conservative_mutation(self) -> bool:
        """Check if residue is functionally equivalent to reference."""
        if self.reference_residue is None:
            return False
        if self.resname in EQUIVALENT_RESIDUES:
            if self.reference_residue.resname in EQUIVALENT_RESIDUES[self.resname]:
                return True
        return False

    @property
    def is_standard(self) -> bool:
        """Check if standard amino acid."""
        return self.resname.upper() in STANDARD_RESIDUES

    @property
    def has_double_funcloc(self) -> bool:
        """Check if residue has multiple function locations."""
        return len(self.funclocs) > 1

    @property
    def has_main_chain_function(self) -> bool:
        """Check if residue has main chain function."""
        for funcloc in self.funclocs:
            if 'main' in funcloc.lower():
                return True
        return False

    @classmethod
    def from_json(
        cls, residue: Dict[str, Any], chiral_id: Optional[int] = None
    ) -> Generator['UniResidue', None, None]:
        """Create UniResidue objects from M-CSA JSON data.

        Args:
            residue: Residue entry from M-CSA JSON
            chiral_id: Stereoisomer identifier. Defaults to None.

        Yields:
            UniResidue objects for each sequence position
        """
        try:
            for uniprot_res in residue['residue_sequences']:
                mcsa_id = residue['mcsa_id']
                uniprot_id = uniprot_res['uniprot_id']
                resname = uniprot_res['code']
                resid = uniprot_res['resid']
                is_reference = uniprot_res['is_reference']
                funclocs = [residue['function_location_abv']]

                yield UniResidue(
                    mcsa_id, uniprot_id, resname, resid, funclocs, is_reference, chiral_id
                )
        except KeyError:
            return


__all__ = ['UniResidue']
