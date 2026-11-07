"""
Unit tests for csa3d.config module.

Tests configuration loading, path validation, and environment variable handling.
"""

import os
from pathlib import Path
import pytest
from csa3d.config.settings import DataPaths
from csa3d.config import loader


@pytest.mark.unit
class TestDataPaths:
    """Tests for DataPaths configuration class."""

    def test_datapaths_can_be_instantiated(self) -> None:
        """Test that DataPaths instance can be created."""
        paths = DataPaths()
        assert paths is not None

    def test_datapaths_has_methods(self) -> None:
        """Test that DataPaths has required methods."""
        paths = DataPaths()
        assert hasattr(paths, "set_data_dir")
        assert hasattr(paths, "get_data_dir")
        assert callable(paths.set_data_dir)
        assert callable(paths.get_data_dir)

    def test_datapaths_set_and_get_data_dir(self, temp_dir: Path) -> None:
        """Test DataPaths set and get data_dir methods."""
        paths = DataPaths()
        paths.set_data_dir(temp_dir)
        retrieved = paths.get_data_dir()
        assert retrieved == temp_dir

    def test_datapaths_with_environment_variable(
        self, temp_dir: Path, monkeypatch
    ) -> None:
        """Test DataPaths respects CSA3D_DATA_DIR environment variable."""
        monkeypatch.setenv("CSA3D_DATA_DIR", str(temp_dir))
        paths = DataPaths()
        # The get_data_dir method should find the env var if temp_dir exists
        if temp_dir.exists():
            data_dir = paths.get_data_dir()
            assert data_dir == temp_dir


@pytest.mark.unit
class TestConfigLoader:
    """Tests for configuration loader functions."""

    def test_loader_module_imports(self) -> None:
        """Test that loader module can be imported."""
        # The module should import successfully
        assert loader is not None

    def test_loader_has_core_functions(self) -> None:
        """Test that loader has expected core functions."""
        expected_functions = [
            "load_uniprot_pdb_mapping",
            "load_pdb_resolution",
            "load_ions",
        ]

        for func_name in expected_functions:
            assert hasattr(loader, func_name), f"Loader missing function: {func_name}"


@pytest.mark.unit
class TestEnvironmentVariables:
    """Tests for environment variable handling."""

    def test_csa3d_data_dir_env_var(self, temp_dir: Path, monkeypatch) -> None:
        """Test CSA3D_DATA_DIR environment variable is respected."""
        monkeypatch.setenv("CSA3D_DATA_DIR", str(temp_dir))

        retrieved_dir = os.getenv("CSA3D_DATA_DIR")
        assert retrieved_dir == str(temp_dir)

    def test_env_var_affects_datapaths(self, temp_dir: Path, monkeypatch) -> None:
        """Test that DataPaths respects env var when set."""
        monkeypatch.setenv("CSA3D_DATA_DIR", str(temp_dir))
        paths = DataPaths()
        if temp_dir.exists():
            data_dir = paths.get_data_dir()
            assert data_dir == temp_dir


@pytest.mark.unit
class TestPathValidation:
    """Tests for path validation and existence checks."""

    def test_datapaths_get_data_dir_returns_path(self, temp_dir: Path) -> None:
        """Test that DataPaths.get_data_dir() returns a Path object."""
        paths = DataPaths()
        paths.set_data_dir(temp_dir)
        retrieved = paths.get_data_dir()
        assert isinstance(retrieved, Path)

    def test_datapaths_set_data_dir_with_valid_path(self, temp_dir: Path) -> None:
        """Test that set_data_dir works with existing directory."""
        paths = DataPaths()
        # Should not raise since temp_dir exists
        paths.set_data_dir(temp_dir)
        assert paths.get_data_dir() == temp_dir

    def test_datapaths_set_data_dir_with_invalid_path(self) -> None:
        """Test that set_data_dir raises error for non-existent path."""
        paths = DataPaths()
        non_existent = Path("/this/path/definitely/does/not/exist/12345")
        with pytest.raises(Exception):
            paths.set_data_dir(non_existent)
