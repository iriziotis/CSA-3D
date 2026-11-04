Installation
=============

CSA-3D requires Python 3.8 or later.

PyPI Installation
-----------------

The easiest way to install CSA-3D::

    pip install csa3d

Installation from Source
------------------------

If you want the development version::

    git clone https://github.com/iriziotis/CSA-3D.git
    cd CSA-3D
    pip install -e .

Installing Optional Dependencies
---------------------------------

For clustering and visualization::

    pip install -e ".[clustering]"

For development (testing, linting, type checking)::

    pip install -e ".[dev]"

For examples and notebooks::

    pip install -e ".[examples]"

Verification
------------

To verify the installation works::

    python -c "import csa3d; print(csa3d.__version__)"

Dependencies
------------

CSA-3D depends on:

- **biopython** (≥1.81) — biological sequence and structure handling
- **numpy** (≥1.21.0) — numerical computing
- **pandas** (≥1.3.0) — data analysis and tables
- **scipy** (≥1.7.0) — scientific computing
- **scikit-learn** (≥0.24.0) — machine learning (clustering)
- **rdkit** (≥2021.09.0) — cheminformatics (ligand similarity)
- **natsort** (≥8.0.0) — natural sorting
- **rmsd** (≥1.4.0) — RMSD calculations

For complete list, see ``pyproject.toml``.
