"""
Pytest configuration and fixtures for CSA-3D tests.

Provides:
- Configuration fixtures (paths, settings)
- Temporary directory fixtures
- Test data fixtures
"""

import os
import sys
import tempfile
from pathlib import Path
from typing import Generator, Optional
import pytest


# Add parent directory to path so we can import modules
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from csa3d.config.settings import DataPaths
except ImportError:
    # Package not fully installed - provide a simple mock
    class DataPaths:
        """Mock DataPaths for testing when package is not installed."""
        def __init__(self):
            self._data_dir = None

        def set_data_dir(self, data_dir: Path) -> None:
            self._data_dir = Path(data_dir)

        @property
        def base_dir(self) -> Path:
            if self._data_dir:
                return self._data_dir
            return Path.home() / ".csa3d" / "data"

        @property
        def pdb_dir(self) -> Path:
            return self.base_dir / "pdb"

        @property
        def uniprot_dir(self) -> Path:
            return self.base_dir / "uniprot"


@pytest.fixture(scope="session")
def csa3d_data_dir() -> Path:
    """
    Provides the CSA-3D data directory path.

    Uses CSA3D_DATA_DIR environment variable or falls back to default location.
    """
    data_dir = os.getenv("CSA3D_DATA_DIR")
    if data_dir:
        return Path(data_dir)

    # Try common locations
    home = Path.home()
    for path in [
        home / "csa3d_data",
        home / ".csa3d" / "data",
        Path("/usr/share/csa3d/data"),
        Path("/opt/csa3d/data"),
    ]:
        if path.exists():
            return path

    # Return the default location (may not exist, which is OK for some tests)
    return home / ".csa3d" / "data"


@pytest.fixture(scope="session")
def data_paths(csa3d_data_dir: Path) -> DataPaths:
    """
    Provides DataPaths configuration for tests.

    Returns:
        DataPaths instance with test data directory
    """
    paths = DataPaths()
    if csa3d_data_dir.exists():
        paths.set_data_dir(csa3d_data_dir)
    return paths


@pytest.fixture(scope="function")
def temp_dir() -> Generator[Path, None, None]:
    """
    Provides a temporary directory that is cleaned up after each test.

    Yields:
        Path to temporary directory
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture(scope="function")
def env_with_data_dir(csa3d_data_dir: Path, monkeypatch) -> None:
    """
    Sets up environment with CSA3D_DATA_DIR for tests.

    Args:
        csa3d_data_dir: Data directory path
        monkeypatch: Pytest monkeypatch fixture
    """
    monkeypatch.setenv("CSA3D_DATA_DIR", str(csa3d_data_dir))


@pytest.fixture(scope="session")
def mcsa_db() -> Optional:
    """
    Provides an initialized Mcsa database instance.

    Note: This requires the CSA-3D data files to be present.
    If data is not available, returns None.

    Returns:
        Mcsa instance or None if data unavailable
    """
    try:
        from csa3d import Mcsa
        db = Mcsa()
        return db
    except Exception as e:
        # Data not available - return None
        # Tests should skip if database is not available
        pytest.skip(f"CSA-3D database not available: {e}")
        return None


@pytest.fixture(scope="function")
def minimal_entry_data() -> dict:
    """
    Provides minimal entry data for testing.

    This is a mock entry structure useful for testing without requiring
    the full M-CSA database.

    Returns:
        Dictionary with minimal entry data structure
    """
    return {
        "mcsa_id": 1,
        "e_number": "1.1.1.1",
        "description": "Alcohol dehydrogenase",
        "pdb_refs": ["1HET"],
        "uniprot_refs": ["P00330"],
        "pdb_sites": [],
        "uni_sites": [],
    }


@pytest.fixture(scope="function")
def sample_pdb_id() -> str:
    """
    Provides a sample PDB ID for testing.

    Returns:
        PDB ID string
    """
    return "1HET"


@pytest.fixture(scope="function")
def sample_uniprot_id() -> str:
    """
    Provides a sample UniProt ID for testing.

    Returns:
        UniProt ID string
    """
    return "P00330"


@pytest.fixture(scope="function")
def sample_mcsa_id() -> int:
    """
    Provides a sample M-CSA entry ID for testing.

    Returns:
        M-CSA entry ID
    """
    return 1


@pytest.fixture(scope="function")
def config_file(temp_dir: Path) -> Path:
    """
    Creates a temporary configuration file for testing.

    Args:
        temp_dir: Temporary directory from fixture

    Returns:
        Path to created config file
    """
    config_path = temp_dir / "config.txt"
    config_path.write_text(
        "[paths]\n"
        "pdb_dir = /tmp/pdb\n"
        "uniprot_file = /tmp/uniprot.fasta\n"
    )
    return config_path


# Markers for conditional test running
def pytest_configure(config):
    """Register custom pytest markers."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "requires_data: marks tests that require CSA-3D data files"
    )
    config.addinivalue_line(
        "markers", "integration: marks tests as integration tests"
    )
    config.addinivalue_line(
        "markers", "unit: marks tests as unit tests"
    )


@pytest.fixture(autouse=True)
def reset_logging():
    """Reset logging configuration between tests."""
    import logging

    # Clear all handlers
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    yield

    # Clean up after test
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)


# Skip tests that require data if data is not available
def pytest_collection_modifyitems(config, items):
    """Automatically skip tests marked with requires_data if data is missing."""
    skip_requires_data = pytest.mark.skip(reason="CSA-3D data not available")

    for item in items:
        if "requires_data" in item.keywords:
            try:
                # Check if data directory exists
                data_dir = os.getenv("CSA3D_DATA_DIR")
                if not data_dir:
                    home = Path.home()
                    data_dir = str(home / ".csa3d" / "data")

                if not Path(data_dir).exists():
                    item.add_marker(skip_requires_data)
            except Exception:
                item.add_marker(skip_requires_data)
