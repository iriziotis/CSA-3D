"""CSA-3D: Extract and analyze enzyme catalytic sites from M-CSA."""

from .data_setup import setup_data
from .config import set_data_dir, get_data_dir

__version__ = "0.2.2"

# Lazy import of main classes to avoid issues with hardcoded paths in root-level config
# These will be imported on first access after setup_data() is called
def __getattr__(name):
    """Lazy import of root-level modules."""
    if name in ("Mcsa", "Entry", "PdbSite", "PdbResidue", "UniResidue", "UniSite", "Superimposer"):
        try:
            import sys
            from pathlib import Path

            # Add parent directory to path for root-level module imports
            _root_dir = str(Path(__file__).parent.parent)
            if _root_dir not in sys.path:
                sys.path.insert(0, _root_dir)

            # Import the requested module
            if name == "Mcsa":
                from Mcsa import Mcsa
                return Mcsa
            elif name == "Entry":
                from Entry import Entry
                return Entry
            elif name == "PdbSite":
                from PdbSite import PdbSite
                return PdbSite
            elif name == "PdbResidue":
                from PdbResidue import PdbResidue
                return PdbResidue
            elif name == "UniResidue":
                from UniResidue import UniResidue
                return UniResidue
            elif name == "UniSite":
                from UniSite import UniSite
                return UniSite
            elif name == "Superimposer":
                from Superimposer import Superimposer
                return Superimposer
        except ImportError as e:
            raise AttributeError(
                f"Could not import {name}. Make sure to run setup_data() first:\n"
                f"  from csa3d import setup_data\n"
                f"  setup_data('/path/to/data')\n"
                f"Original error: {e}"
            )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "setup_data",
    "set_data_dir",
    "get_data_dir",
    "Mcsa",
    "Entry",
    "PdbSite",
    "PdbResidue",
    "UniResidue",
    "UniSite",
    "Superimposer",
]
