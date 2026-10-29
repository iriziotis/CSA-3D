"""Data loading utilities for CSA-3D.

Loads and caches data files required by CSA-3D.
"""

import csv
import json
import warnings
from collections import defaultdict
from glob import glob
from pathlib import Path
from typing import Any, Dict, Set, Tuple

from .settings import _data_paths


def _load_json_mapping(file_path: Path) -> Dict[Any, Any]:
    """Load a JSON file as a dictionary.

    Args:
        file_path: Path to JSON file

    Returns:
        Loaded JSON data

    Raises:
        FileNotFoundError: If file does not exist
        ValueError: If JSON is invalid
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Data file not found: {file_path}")

    try:
        with open(file_path, "r") as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in {file_path}: {e}")


def load_uniprot_pdb_mapping() -> Dict[str, list]:
    """Load UniProt to PDB mapping from SIFTS.

    Returns:
        Dictionary mapping UniProt IDs to list of PDB IDs

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("uniprot_pdb_mapping")

    if not path.exists():
        raise FileNotFoundError(
            f"UniProt-PDB mapping file not found: {path}\n"
            "Run csa3d.config.setup_data() to download"
        )

    mapping = {}
    with open(path, "r") as f:
        next(f)  # Skip header
        for line in csv.reader(f):
            if len(line) >= 2:
                mapping[line[0]] = line[1].split(";")

    return mapping


def load_pdb_uniprot_mapping() -> Dict[Tuple[str, str], str]:
    """Load PDB to UniProt mapping from SIFTS.

    Returns:
        Dictionary mapping (PDB ID, chain) tuples to UniProt IDs

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("pdb_uniprot_mapping")

    if not path.exists():
        raise FileNotFoundError(
            f"PDB-UniProt mapping file not found: {path}\n"
            "Run csa3d.config.setup_data() to download"
        )

    mapping = {}
    with open(path, "r") as f:
        next(f)  # Skip header
        next(f)  # Skip second header line
        for line in csv.reader(f):
            if len(line) >= 3:
                mapping[(line[0], line[1])] = line[2]

    return mapping


def load_pdb_ec_mapping() -> Dict[Tuple[str, str], str]:
    """Load PDB to EC number mapping from SIFTS.

    Returns:
        Dictionary mapping (PDB ID, chain) tuples to EC numbers

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("pdb_ec_mapping")

    if not path.exists():
        raise FileNotFoundError(f"PDB-EC mapping file not found: {path}")

    mapping = {}
    with open(path, "r") as f:
        next(f)  # Skip header
        for line in csv.reader(f):
            if len(line) >= 4:
                mapping[(line[0], line[1])] = line[3]

    return mapping


def load_kegg_reactions() -> Dict[str, list]:
    """Load KEGG reaction components mapping.

    Returns:
        Dictionary mapping EC numbers to reaction components

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("kegg_reactions")

    if not path.exists():
        raise FileNotFoundError(f"KEGG reactions file not found: {path}")

    mapping = {}
    with open(path, "r") as f:
        for line in csv.reader(f):
            if len(line) >= 2:
                mapping[line[0]] = [line[1].split(";")]

    return mapping


def load_rhea_reactions() -> Dict[str, list]:
    """Load RHEA reaction components mapping.

    Returns:
        Dictionary mapping EC numbers to reaction components

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("rhea_reactions")

    if not path.exists():
        raise FileNotFoundError(f"RHEA reactions file not found: {path}")

    mapping = {}
    with open(path, "r") as f:
        next(f)  # Skip header
        for line in csv.reader(f):
            if len(line) >= 4:
                mapping[line[2]] = line[3].split(";")

    return mapping


def load_crystallization_artifacts() -> Set[str]:
    """Load set of crystallization artifact HET codes.

    Returns:
        Set of HET codes that are crystallization artifacts

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("crystallization_hets")

    if not path.exists():
        raise FileNotFoundError(f"Crystallization HETs file not found: {path}")

    artifacts = set()
    with open(path, "r") as f:
        for line in f:
            artifacts.add(line.strip())

    return artifacts


def load_ions() -> Tuple[Set[str], Set[str], Set[str]]:
    """Load ion classifications.

    Returns:
        Tuple of (metals, reactive_nonmetals, noble_gases) sets

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("pdb_ions")

    if not path.exists():
        raise FileNotFoundError(f"PDB ions file not found: {path}")

    metals = set()
    reactive_nonmetals = set()
    noble_gases = set()

    with open(path, "r") as f:
        next(f)  # Skip header
        for line in csv.reader(f):
            if len(line) >= 4:
                ion_type = line[3]
                if ion_type == "METAL":
                    metals.add(line[0])
                elif ion_type == "REACTIVE_NONMETAL":
                    reactive_nonmetals.add(line[0])
                elif ion_type == "NOBLE_GAS":
                    noble_gases.add(line[0])

    return metals, reactive_nonmetals, noble_gases


def load_pdb_cofactors() -> Dict[str, Set[str]]:
    """Load PDB ID to cofactors mapping.

    Returns:
        Dictionary mapping PDB IDs to sets of cofactors

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("pdb_cofactors")

    if not path.exists():
        raise FileNotFoundError(f"PDB cofactors file not found: {path}")

    cofactors = defaultdict(set)
    with open(path, "r") as f:
        next(f)  # Skip header
        for line in csv.reader(f, quotechar='"'):
            if len(line) >= 5:
                cofactors[line[0]].add(line[4])

    return dict(cofactors)


def load_mcsa_cofactors() -> Dict[str, Set[str]]:
    """Load M-CSA EC number to cofactors mapping.

    Returns:
        Dictionary mapping EC numbers to sets of cofactors

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("mcsa_cofactors")

    if not path.exists():
        raise FileNotFoundError(f"M-CSA cofactors file not found: {path}")

    cofactors = defaultdict(set)
    with open(path, "r") as f:
        next(f)  # Skip header
        for line in csv.reader(f, quotechar='"'):
            if len(line) >= 5:
                for ec in line[4].split(";"):
                    cofactors[ec].add(line[1])

    return dict(cofactors)


def load_pdb_resolution() -> Dict[str, str]:
    """Load PDB resolution data.

    Returns:
        Dictionary mapping PDB IDs to resolution values

    Raises:
        FileNotFoundError: If data file not found
    """
    path = _data_paths.get_path("pdb_resolution")

    if not path.exists():
        raise FileNotFoundError(f"PDB resolution file not found: {path}")

    resolution = {}
    with open(path, "r") as f:
        for line in csv.reader(f, delimiter=";"):
            if len(line) >= 2:
                resolution[line[0].strip().lower()] = line[1].strip()

    return resolution
