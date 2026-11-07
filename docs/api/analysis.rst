Structural Analysis
===================

Superposition and comparison of catalytic sites.

Superimposer
------------

.. autoclass:: csa3d.Superimposer
   :members:
   :undoc-members:
   :show-inheritance:

The Superimposer class implements the hybrid iterative and weighted superposition
algorithm described in the CSA-3D paper. It handles structural alignment accounting
for conformational variation in catalytic sites.

Key Methods
^^^^^^^^^^^

- :meth:`~csa3d.Superimposer.set`: Configure coordinates and fitting parameters
- :meth:`~csa3d.Superimposer.run_unweighted`: Fast iterative fitting with outlier rejection
- :meth:`~csa3d.Superimposer.run_weighted`: Precise Gaussian-weighted fitting
- :meth:`~csa3d.Superimposer.get_rotran`: Retrieve rotation and translation
- :meth:`~csa3d.Superimposer.get_rms`: Retrieve RMSD values

Example
^^^^^^^

.. code-block:: python

    from csa3d import Superimposer
    import numpy as np
    
    # Create coordinate sets
    sup = Superimposer()
    ref_coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    mob_coords = np.array([[0.1, 0, 0], [1.1, 0.1, 0], [0, 1.1, 0]], dtype=float)
    
    # Configure and fit
    sup.set(ref_coords, mob_coords, cycles=3, cutoff=3.0)
    sup.run_weighted()
    
    # Get results
    rot, tran = sup.get_rotran()
    rms, rms_all = sup.get_rms()
    print(f"RMSD: {rms:.2f} Å")
