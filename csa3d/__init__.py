"""CSA-3D: Extract and analyze enzyme catalytic sites from M-CSA."""

from .data_setup import setup_data
from .config import set_data_dir, get_data_dir

__version__ = "0.2.2"

# Add parent directory to path for root-level module imports
import sys
from pathlib import Path
_root_dir = str(Path(__file__).parent.parent)
if _root_dir not in sys.path:
    sys.path.insert(0, _root_dir)

# Main classes will be imported after setup_data() and set_data_dir() are called
# This is declared in __all__ for documentation purposes
__all__ = [
    "setup_data",
    "set_data_dir",
    "get_data_dir",
]

def _import_main_classes():
    """Import main classes after data is configured."""
    global Mcsa, Entry, PdbSite, PdbResidue, UniResidue, UniSite, Superimposer
    try:
        from Mcsa import Mcsa
        from Entry import Entry
        from PdbSite import PdbSite
        from PdbResidue import PdbResidue
        from UniResidue import UniResidue
        from UniSite import UniSite
        from Superimposer import Superimposer

        # Add to module globals so they're accessible
        globals()['Mcsa'] = Mcsa
        globals()['Entry'] = Entry
        globals()['PdbSite'] = PdbSite
        globals()['PdbResidue'] = PdbResidue
        globals()['UniResidue'] = UniResidue
        globals()['UniSite'] = UniSite
        globals()['Superimposer'] = Superimposer

        __all__.extend([
            "Mcsa", "Entry", "PdbSite", "PdbResidue",
            "UniResidue", "UniSite", "Superimposer"
        ])
        return True
    except Exception as e:
        return False

# Try to import on module load
_import_main_classes()

def __getattr__(name):
    """Provide helpful error message if classes haven't been imported yet."""
    if name in ("Mcsa", "Entry", "PdbSite", "PdbResidue", "UniResidue", "UniSite", "Superimposer"):
        raise AttributeError(
            f"\nCould not import {name}.\n\n"
            f"Make sure to set up data first:\n"
            f"  from csa3d import setup_data, set_data_dir\n"
            f"  setup_data('/path/to/data')\n"
            f"  set_data_dir('/path/to/data')\n\n"
            f"Then restart your Python session and try again.\n"
        )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
