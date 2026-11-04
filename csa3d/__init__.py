"""CSA-3D: extract and analyse enzyme catalytic sites in 3D from the M-CSA.

    >>> import csa3d
    >>> csa3d.setup_data()          # one-off download to ~/.csa3d/data
    >>> db = csa3d.Mcsa()
    >>> db.build(entry_ids=1)
"""

from .config import get_data_dir, set_data_dir, validate_config
from .data_setup import setup_data
from .Entry import Entry
from .Mcsa import Mcsa
from .PdbResidue import Het, PdbResidue
from .PdbSite import PdbSite
from .Superimposer import Superimposer
from .UniResidue import UniResidue
from .UniSite import UniSite

__version__ = "0.4.1"

__all__ = [
    "setup_data",
    "set_data_dir",
    "get_data_dir",
    "validate_config",
    "Mcsa",
    "Entry",
    "PdbSite",
    "PdbResidue",
    "Het",
    "UniSite",
    "UniResidue",
    "Superimposer",
]
