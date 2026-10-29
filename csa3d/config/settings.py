"""Configuration settings for CSA-3D with environment variable support."""

import os
from pathlib import Path
from typing import Dict, Optional

# Default data directory locations to try
_DEFAULT_PATHS = [
    "/Users/riziotis/ebi/phd/data",
    "/nfs/research/thornton/riziotis/research/phd/data",
    str(Path.home() / ".csa3d" / "data"),
]


class ConfigError(Exception):
    """Raised when configuration is invalid or incomplete."""

    pass


class DataPaths:
    """Manages data file paths for CSA-3D.

    Paths can be configured via:
    1. Environment variable CSA3D_DATA_DIR
    2. Direct call to set_data_dir()
    3. Default locations
    """

    def __init__(self) -> None:
        """Initialize data paths."""
        self._data_dir: Optional[Path] = None
        self._initialized = False

    def set_data_dir(self, data_dir: str) -> None:
        """Set the data directory explicitly.

        Args:
            data_dir: Path to directory containing CSA-3D data files

        Raises:
            ConfigError: If directory does not exist
        """
        path = Path(data_dir)
        if not path.is_dir():
            raise ConfigError(f"Data directory does not exist: {data_dir}")
        self._data_dir = path
        self._initialized = True

    def get_data_dir(self) -> Path:
        """Get the configured data directory.

        Returns:
            Path to data directory

        Raises:
            ConfigError: If no valid data directory is configured
        """
        if self._data_dir is not None:
            return self._data_dir

        # Try environment variable first
        env_dir = os.environ.get("CSA3D_DATA_DIR")
        if env_dir:
            path = Path(env_dir)
            if path.is_dir():
                self._data_dir = path
                self._initialized = True
                return path
            raise ConfigError(
                f"CSA3D_DATA_DIR points to non-existent directory: {env_dir}"
            )

        # Try default locations
        for default_path in _DEFAULT_PATHS:
            path = Path(default_path)
            if path.is_dir():
                self._data_dir = path
                self._initialized = True
                return path

        # No valid path found
        raise ConfigError(
            "No CSA-3D data directory found. "
            "Please either:\n"
            "  1. Set environment variable: export CSA3D_DATA_DIR=/path/to/data\n"
            "  2. Call: csa3d.config.set_data_dir('/path/to/data')\n"
            "  3. Run: csa3d.config.setup_data() to download datasets\n"
            f"Searched in: {', '.join(_DEFAULT_PATHS)}"
        )

    def validate(self) -> bool:
        """Validate that all required data files exist.

        Returns:
            True if configuration is valid

        Raises:
            ConfigError: If required files are missing
        """
        data_dir = self.get_data_dir()
        required_files = [
            "mcsa/catalytic_residues_homologues.json",
            "mcsa/catalytic_residues_references.json",
            "sifts/uniprot_pdb.csv",
            "sifts/pdb_chain_uniprot.csv",
            "sifts/pdb_chain_enzyme.csv",
        ]

        missing = []
        for rel_path in required_files:
            full_path = data_dir / rel_path
            if not full_path.exists():
                missing.append(rel_path)

        if missing:
            raise ConfigError(
                f"Missing required data files in {data_dir}:\n"
                + "\n".join(f"  - {f}" for f in missing)
                + "\nRun: csa3d.config.setup_data() to download"
            )
        return True

    def get_path(self, key: str) -> Path:
        """Get path to a specific data file.

        Args:
            key: File identifier (e.g., 'catalytic_residues', 'pdb_assemblies')

        Returns:
            Path to the file

        Raises:
            ConfigError: If file key is unknown
        """
        data_dir = self.get_data_dir()

        paths = {
            "catalytic_residues": data_dir / "mcsa" / "catalytic_residues_homologues.json",
            "catalytic_residues_references": data_dir / "mcsa" / "catalytic_residues_references.json",
            "mcsa_entry_info": data_dir / "mcsa" / "entry_info" / "mcsa_entry_info.*.json",
            "pdb_assemblies": data_dir / "pdbe" / "assembly_cif" / "cif",
            "kegg_molecules": data_dir / "kegg" / "sdf",
            "chebi_molecules": data_dir / "chebi" / "sdf",
            "pdbe_het_molecules": data_dir / "pdbe" / "het_sdf",
            "uniprot_pdb_mapping": data_dir / "sifts" / "uniprot_pdb.csv",
            "pdb_uniprot_mapping": data_dir / "sifts" / "pdb_chain_uniprot.csv",
            "pdb_ec_mapping": data_dir / "sifts" / "pdb_chain_enzyme.csv",
            "kegg_reactions": data_dir / "kegg" / "kegg_reaction_components.csv",
            "rhea_reactions": data_dir / "rhea" / "rhea_reaction_components.csv",
            "pdb_ions": data_dir / "pdbe" / "pdb_ions.csv",
            "pdb_cofactors": data_dir / "pdbe" / "pdb_cofactors.csv",
            "mcsa_cofactors": data_dir / "mcsa" / "mcsa_cofactors.csv",
            "crystallization_hets": data_dir / "pdbe" / "crystallization_hets.csv",
            "pdb_resolution": data_dir / "pdbe" / "resolu.idx",
        }

        if key not in paths:
            raise ConfigError(
                f"Unknown path key: {key}. "
                f"Valid keys: {', '.join(paths.keys())}"
            )

        return paths[key]


# Global instance
_data_paths = DataPaths()


def get_config() -> Dict[str, Path]:
    """Get current configuration as dictionary.

    Returns:
        Dictionary of configuration values
    """
    data_dir = _data_paths.get_data_dir()
    return {
        "data_dir": data_dir,
        "catalytic_residues": _data_paths.get_path("catalytic_residues"),
        "pdb_assemblies": _data_paths.get_path("pdb_assemblies"),
    }


def set_data_dir(data_dir: str) -> None:
    """Set the CSA-3D data directory.

    Args:
        data_dir: Path to directory containing CSA-3D data files

    Raises:
        ConfigError: If directory does not exist
    """
    _data_paths.set_data_dir(data_dir)


def get_data_dir() -> str:
    """Get the current CSA-3D data directory.

    Returns:
        Path to data directory as string

    Raises:
        ConfigError: If no valid configuration found
    """
    return str(_data_paths.get_data_dir())


def validate_config() -> bool:
    """Validate current configuration.

    Returns:
        True if all required data files exist

    Raises:
        ConfigError: If configuration is invalid
    """
    return _data_paths.validate()
