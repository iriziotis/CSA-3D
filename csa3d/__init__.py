"""CSA-3D: Extract and analyze enzyme catalytic sites from M-CSA."""

from .data_setup import setup_data
from .config import set_data_dir, get_data_dir

__version__ = "0.2.2"

__all__ = [
    "setup_data",
    "set_data_dir",
    "get_data_dir",
]

# Note: The Mcsa class is available at the package root level.
# Import with: from Mcsa import Mcsa (after setting up data with csa3d.setup_data())
