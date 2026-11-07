"""Parsers for the CSA-3D data files.

Each function reads one file from the configured data directory. Use
``csa3d.datasets`` for cached access; these functions re-read on every call.
"""

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Set, Tuple

from .settings import _data_paths


def _require(key: str) -> Path:
    path = _data_paths.get_path(key)
    if not path.exists():
        raise FileNotFoundError(
            f"Data file not found: {path}\n"
            "Run csa3d.setup_data() to download the CSA-3D data, or point "
            "csa3d.set_data_dir() to an existing data directory."
        )
    return path


def _rows(path: Path, skip_header: bool = True, **reader_kwargs) -> Iterator[List[str]]:
    """Yield CSV rows, skipping '#' comment lines (SIFTS) and the header row."""
    with open(path, "r") as f:
        lines = (line for line in f if not line.startswith("#"))
        reader = csv.reader(lines, **reader_kwargs)
        if skip_header:
            next(reader, None)
        yield from reader


def load_uniprot_pdb_mapping() -> Dict[str, List[str]]:
    """UniProt accession -> list of PDB IDs (SIFTS)."""
    return {
        row[0]: row[1].split(";")
        for row in _rows(_require("uniprot_pdb_mapping"))
        if len(row) >= 2
    }


def load_pdb_uniprot_mapping() -> Dict[Tuple[str, str], str]:
    """(PDB ID, chain) -> UniProt accession (SIFTS)."""
    return {
        (row[0], row[1]): row[2]
        for row in _rows(_require("pdb_uniprot_mapping"))
        if len(row) >= 3
    }


def load_pdb_ec_mapping() -> Dict[Tuple[str, str], str]:
    """(PDB ID, chain) -> EC number (SIFTS)."""
    return {
        (row[0], row[1]): row[3]
        for row in _rows(_require("pdb_ec_mapping"))
        if len(row) >= 4
    }


def load_kegg_reactions() -> Dict[str, List[str]]:
    """EC number -> KEGG compound IDs."""
    return {
        row[0]: row[1].split(";")
        for row in _rows(_require("kegg_reactions"))
        if len(row) >= 2
    }


def load_rhea_reactions() -> Dict[str, List[str]]:
    """EC number -> ChEBI IDs of one RHEA reaction for that EC."""
    return {
        row[2]: row[3].split(";")
        for row in _rows(_require("rhea_reactions"))
        if len(row) >= 4
    }


def load_crystallization_artifacts() -> Set[str]:
    """HET codes of common crystallization artefacts."""
    with open(_require("crystallization_hets"), "r") as f:
        return {line.strip() for line in f if line.strip()}


def load_ions() -> Tuple[Set[str], Set[str], Set[str]]:
    """(metals, reactive non-metals, noble gases) HET code sets."""
    metals: Set[str] = set()
    reactive_nonmetals: Set[str] = set()
    noble_gases: Set[str] = set()
    groups = {
        "METAL": metals,
        "REACTIVE_NONMETAL": reactive_nonmetals,
        "NOBLE_GAS": noble_gases,
    }
    for row in _rows(_require("pdb_ions")):
        if len(row) >= 4 and row[3] in groups:
            groups[row[3]].add(row[0])
    return metals, reactive_nonmetals, noble_gases


def load_pdb_cofactors() -> Dict[str, Set[str]]:
    """PDB ID -> HET codes of cofactors bound in that entry."""
    cofactors: Dict[str, Set[str]] = defaultdict(set)
    for row in _rows(_require("pdb_cofactors"), quotechar='"'):
        if len(row) >= 5:
            cofactors[row[0]].add(row[4])
    return cofactors


def load_mcsa_cofactors() -> Dict[str, Set[str]]:
    """EC number -> HET codes of cofactors annotated in M-CSA."""
    cofactors: Dict[str, Set[str]] = defaultdict(set)
    for row in _rows(_require("mcsa_cofactors"), quotechar='"'):
        if len(row) >= 5:
            for ec in row[4].split(";"):
                cofactors[ec].add(row[1])
    return cofactors


def load_pdb_resolution() -> Dict[str, str]:
    """Lower-case PDB ID -> resolution string (wwPDB resolu.idx)."""
    return {
        row[0].strip().lower(): row[1].strip()
        for row in _rows(_require("pdb_resolution"), skip_header=False, delimiter=";")
        if len(row) >= 2
    }
