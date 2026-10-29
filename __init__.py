"""CSA-3D: Extract and analyze enzyme catalytic sites in 3D.

This package provides tools to extract, analyze, and manipulate enzyme
catalytic sites from the Mechanism and Catalytic Site Atlas (M-CSA).

Main usage:
    from csa3d import Mcsa
    db = Mcsa()
    db.build([1, 2, 3])  # Build entries 1, 2, 3
    for entry in db:
        print(entry.mcsa_id, len(entry.pdbsites))

See documentation for more details: https://csa3d.readthedocs.io
"""

from .Mcsa import Mcsa
from csa3d_version import __version__

__all__ = ["Mcsa", "__version__"]
