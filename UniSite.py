"""UniProt catalytic site representation for M-CSA."""

from typing import Optional, List, Dict, Any, Generator

# Handle imports from both package and root-level contexts
try:
    from residue_definitions import AA_3TO1
    from UniResidue import UniResidue
except ImportError:
    from .residue_definitions import AA_3TO1
    from .UniResidue import UniResidue


class UniSite:
    """Catalytic site in UniProt sequence coordinates.

    Represents a catalytic site at the sequence level with residue positions
    and mappings to PDB structures.

    Attributes:
        parent_entry: Parent Entry object
        residues: List of UniResidue objects
        residues_dict: Dictionary mapping full_id to UniResidue
        mapped_pdbsites: List of mapped PdbSite objects
        reference_site: Reference site for alignment
    """

    def __init__(self) -> None:
        """Initialize UniSite."""
        self.parent_entry: Optional[Any] = None
        self.residues: List[UniResidue] = []
        self.residues_dict: Dict[tuple, UniResidue] = {}
        self.mapped_pdbsites: List[Any] = []
        self.reference_site: Optional['UniSite'] = None

    def __str__(self) -> str:
        """Return pseudo-sequence in one-letter code."""
        return self.sequence

    def __iter__(self) -> Generator[UniResidue, None, None]:
        """Iterate over residues."""
        yield from self.residues

    def __contains__(self, residue: UniResidue) -> bool:
        """Check if residue is in this site."""
        return residue.full_id in self.residues_dict

    def __getitem__(self, full_id: tuple) -> UniResidue:
        """Return residue with given full_id."""
        return self.residues_dict[full_id]

    # Alternative constructors

    @classmethod
    def build_reference(cls, reslist: List[UniResidue]) -> 'UniSite':
        """Build reference site from UniResidue list.

        Args:
            reslist: List of reference UniResidue objects

        Returns:
            Reference UniSite object
        """
        ref = UniSite.from_list(reslist)
        ref.reference_site = ref
        return ref

    @classmethod
    def build(
        cls, reslist: List[UniResidue], reference_site: Optional['UniSite'] = None
    ) -> 'UniSite':
        """Build UniSite and map to reference.

        Args:
            reslist: List of UniResidue objects
            reference_site: Reference site for alignment. Defaults to None.

        Returns:
            UniSite object
        """
        site = cls.from_list(reslist)
        if site.is_reference:
            site.reference_site = site
        else:
            site.reference_site = reference_site
            site._map_reference_residues()
        return site

    @classmethod
    def from_list(cls, reslist: List[UniResidue]) -> 'UniSite':
        """Construct UniSite from residue list.

        Args:
            reslist: List of UniResidue objects

        Returns:
            UniSite object
        """
        site = cls()
        reslist = UniSite._cleanup_list(reslist)
        for res in reslist:
            site.add(res)
        return site

    # Properties

    @property
    def mcsa_id(self) -> Optional[int]:
        """Get M-CSA ID from first residue."""
        for res in self.residues:
            if res.mcsa_id:
                return res.mcsa_id
        return None

    @property
    def uniprot_id(self) -> Optional[str]:
        """Get UniProt accession."""
        if self.residues:
            return self.residues[0].uniprot_id
        return None

    @property
    def id(self) -> Optional[str]:
        """Get UniProt site ID (same as uniprot_id)."""
        return self.uniprot_id

    @property
    def sequence(self) -> str:
        """Get pseudo-sequence in one-letter code."""
        return ''.join(
            [
                AA_3TO1.get(res.resname, 'X') if (res.is_standard or res.is_gap) else 'X'
                for res in self
            ]
        )

    @property
    def size(self) -> int:
        """Get site size in residue count."""
        return len(self.residues)

    @property
    def acts_on_polymer(self) -> bool:
        """Check if enzyme substrate is polymer (protein/nucleic acid)."""
        if self.parent_entry is None:
            return False
        return self.parent_entry.info['reaction']['is_polymeric']

    @property
    def is_reference(self) -> bool:
        """Check if this is reference site."""
        if self.size > 0:
            return self.residues[0].is_reference
        return False

    @property
    def is_conserved(self) -> bool:
        """Check if all residues conserved vs reference."""
        if self.reference_site is None:
            return self.is_reference
        return str(self) == str(self.reference_site) or self.is_reference

    @property
    def is_conservative_mutation(self) -> bool:
        """Check if mutations are conservative."""
        result = False
        for res in self.residues:
            if not res.is_conserved and not res.is_conservative_mutation:
                return False
            if res.is_conservative_mutation:
                result = True
        return result

    def add(self, residue: UniResidue) -> bool:
        """Add UniResidue to site.

        Args:
            residue: UniResidue to add

        Returns:
            True if successful
        """
        if not isinstance(residue, UniResidue):
            return False
        self.residues.append(residue)
        self.residues_dict[residue.full_id] = residue
        return True

    # Private methods

    def _map_reference_residues(self) -> None:
        """Map residues to reference positions and fill gaps."""
        if self.reference_site is None:
            return
        for reference_residue in self.reference_site:
            found = False
            for res in self:
                if reference_residue == res.reference_residue:
                    found = True
                    break
            if not found:
                gap = UniResidue()
                gap.reference_residue = reference_residue
                self.add(gap)
        self._reorder()

    def _reorder(self) -> None:
        """Reorder residues to match reference ordering."""
        if self.reference_site is None:
            return
        reorder = []
        for i, reference_residue in enumerate(self.reference_site):
            for j, res in enumerate(self):
                if reference_residue == res.reference_residue:
                    reorder.append(j)
                    break
        if len(reorder) == len(self.residues):
            self.residues = [self.residues[i] for i in reorder]

    @staticmethod
    def _cleanup_list(reslist: List[UniResidue]) -> List[UniResidue]:
        """Remove duplicate residues with different funclocs.

        Args:
            reslist: List of UniResidue objects

        Returns:
            Cleaned list without redundant residues
        """
        new_reslist: List[UniResidue] = []
        seen: set = set()
        ignore: set = set()
        for p in reslist:
            for q in reslist:
                if p == q or (q.id, p.id) in seen:
                    continue
                if p.is_equivalent(q, by_chiral_id=False):
                    if p.funclocs != q.funclocs:
                        new_res = p.copy()
                        new_res.funclocs = [p.funclocs[0], q.funclocs[0]]
                        new_reslist.append(new_res)
                        ignore.add(p.id)
                        ignore.add(q.id)
                seen.add((p.id, q.id))
        for p in reslist:
            if p.id not in ignore and p not in new_reslist:
                new_reslist.append(p)
        return new_reslist


__all__ = ['UniSite']
