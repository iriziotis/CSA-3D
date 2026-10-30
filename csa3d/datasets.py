"""Lazily loaded reference datasets used throughout CSA-3D.

Attributes are resolved on first access against the currently configured data
directory (see ``csa3d.set_data_dir``) and cached per directory, so importing
csa3d never touches the data files.

    >>> from csa3d import datasets
    >>> datasets.PDB2EC[("1b73", "A")]
"""

from typing import Any, Callable, Dict, Tuple

from .config import loader
from .config.settings import _data_paths

_PATHS = {
    "ASSEMBLIES_DIR": "pdb_assemblies",
    "KEGG_MOLS_DIR": "kegg_molecules",
    "CHEBI_MOLS_DIR": "chebi_molecules",
    "HET_MOLS_DIR": "pdbe_het_molecules",
    "CAT_RES_INFO": "catalytic_residues",
    "REF_RES_INFO": "catalytic_residues_references",
    "MCSA_ENTRY_INFO": "mcsa_entry_info",
}

_LOADERS: Dict[str, Callable[[], Any]] = {
    "UNI2PDB": loader.load_uniprot_pdb_mapping,
    "PDB2UNI": loader.load_pdb_uniprot_mapping,
    "PDB2EC": loader.load_pdb_ec_mapping,
    "KEGG_EC_REACTION": loader.load_kegg_reactions,
    "RHEA_EC_REACTION": loader.load_rhea_reactions,
    "CRYSTALLIZATION_HETS": loader.load_crystallization_artifacts,
    "METALS": lambda: loader.load_ions()[0],
    "REACTIVE_NONMETALS": lambda: loader.load_ions()[1],
    "NOBLE_GASES": lambda: loader.load_ions()[2],
    "PDBID_COFACTORS": loader.load_pdb_cofactors,
    "MCSA_EC_COFACTORS": loader.load_mcsa_cofactors,
    "PDBID_RESOLUTION": loader.load_pdb_resolution,
}

_cache: Dict[Tuple[str, str], Any] = {}


def __getattr__(name: str) -> Any:
    if name in _PATHS:
        return str(_data_paths.get_path(_PATHS[name]))
    if name in _LOADERS:
        key = (str(_data_paths.get_data_dir()), name)
        if key not in _cache:
            _cache[key] = _LOADERS[name]()
        return _cache[key]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(_PATHS) + list(_LOADERS))
