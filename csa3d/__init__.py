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

# Try to import main classes
# These are available if data has been set up
Mcsa = None
Entry = None
PdbSite = None
PdbResidue = None
UniResidue = None
UniSite = None
Superimposer = None

try:
    from Mcsa import Mcsa
    from Entry import Entry
    from PdbSite import PdbSite
    from PdbResidue import PdbResidue
    from UniResidue import UniResidue
    from UniSite import UniSite
    from Superimposer import Superimposer
except Exception:
    # Classes not available yet - will show error when user tries to use them
    pass

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
