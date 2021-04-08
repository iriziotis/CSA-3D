import numpy as np
from numpy import dot, transpose, sqrt
from numpy.linalg import svd, det, norm
from copy import copy

class Superimposer:

    def __init__(self):
        """Initialize class."""
        self._clear()

    # Public methods

    def set(self, reference_coords, coords, cycles=1, cutoff=999, scaling_factor=None):
        """Set coordinates and parameters for superposition.
            - reference_coords: NxD NumPy array (N: num of points, D: dimensions)
            - coords:  NxD array
            - cycles: Number of iterations of initial unweighted 
                      superposition with outlier rejectioa. If None, no outlier
                      rejection is performed
            - cutoff: Distance cutoff to reject outlying atoms during
                      the initial superposition cycles. If None, no outlier rejection
                      is performed
            - scaling_factor: Arbitrary factor to adjust the Gaussian weighting
                              function. 
        For details, see https://dx.doi.org/10.1529%2Fbiophysj.105.066654
        """
        # Reinitialize
        self._clear()
        # Set coords
        self.reference_coords = reference_coords
        self.coords = coords
        # Check if sets have the same size
        n = reference_coords.shape
        m = coords.shape
        if n != m:
            raise Exception("Coordinate sets do not have the same dimensions.")
        self.n = n[0]
        # Set superposition parameters
        self.cycles = cycles 
        self.cutoff = cutoff
        if scaling_factor:
            self.scaling_factor = scaling_factor
        else:
            self.scaling_factor = self._auto_scaling_factor()


    def run_unweighted(self):
        """Classic iterative superposition. After each iteration, residues
        with an average distance higher than the specified 'cutoff' are rejected. Number
        of iterations is defined in 'cycles'.
        """
        reference_coords = self.reference_coords.copy()
        coords = self.coords.copy()
        coords_all = self.coords.copy()
        untransformed_coords = self.coords.copy()
        min_rms = 999

        # Run superposition
        for i in range(self.cycles):
            # Fit
            rot, tran = self._fit(coords, reference_coords)
            # Transform coords
            coords = self._transform(coords, rot, tran)
            coords_all = self._transform(coords_all, rot, tran)
            # Calculate RMSDs
            rms = self._rms(reference_coords, coords)
            rms_all = self._rms(self.reference_coords, coords_all)
            # Reject outliers
            if rms < min_rms:
                min_rms = rms
            diff = norm(reference_coords-coords, axis=1)
            to_keep = np.where(diff < self.cutoff)
            # We always want at least three atoms to be superimposed
            if to_keep[0].size < 3:
                break
            coords = coords[to_keep]
            reference_coords = reference_coords[to_keep]
            untransformed_coords = untransformed_coords[to_keep]

        # Final rotation matrix, translation vector and RMSD
        self.rot, self.tran = self._fit(untransformed_coords, coords)
        self.rms = min_rms
        self.rms_all = rms_all
        self.temp_coords = coords.copy() 
        self.temp_reference_coords = reference_coords.copy()
        self.temp_untransformed_coords = untransformed_coords.copy()


    def run_weighted(self):
        """Hybrid iterative/weighted superposition. First coordinates are superimposed
        iteratively with outlier rejection (user specified cycles and distance cutoff),
        and then the cycle includes multiple nested iterations of weighted superposition
        until RMSD convergence is achieved. If no cycles or cutoff are specified, only
        the second step is performed (no initial outlier rejection iterations).
        """
        coords_all = self.coords.copy()
        n = self.n
        c = self.scaling_factor

        # Initial superposition
        self.run_unweighted()
        coords = self.temp_coords.copy()
        reference_coords = self.temp_reference_coords.copy()
        untransformed_coords = self.temp_untransformed_coords.copy()
        coords_all = self._transform(coords_all, self.rot, self.tran)

        # Weighted iterative superposition
        prev_rms = self.rms
        for i in range(5000):
            # Get weights
            weights = self._weights(coords, reference_coords, c)
            # Weighted fit
            rot, tran = self._fit(coords, reference_coords, weights)
            # Transform
            coords = self._transform(coords, rot, tran)
            coords_all = self._transform(coords_all, rot, tran)
            # Weighted RMSD 
            rms = self._rms(coords, reference_coords, weights)
            # RMSD over all atoms, unweighted
            rms_all = self._rms(self.reference_coords, coords_all)
            # Check for convergence
            if 0 <= prev_rms - rms <= 0.000001:
                break
            prev_rms = rms

        if i == 4999:
            raise Exception("Could not reach convergence after 5000 cycles. Maybe try unweighted fit.")
        else:
            # Final rotation matrix, translation vector and RMSD
            self.rot, self.tran = self._fit(untransformed_coords, coords)
            self.rms = rms
            self.rms_all = rms_all
        return
    
    def get_rotran(self):
        """Right multiplying rotation matrix and translation."""
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rot, self.tran

    def get_rms(self):
        """Root mean square deviation of superimposed coordinates."""
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rms, self.rms_all

    # Private methods

    def _clear(self):
        self.reference_coords = None
        self.coords = None
        self.rot = None
        self.tran = None
        self.rms = None
        self.cycles = None
        self.cutoff = None
        self.scaling_factor = None

    def _rms(self, p_coords, q_coords, weights=None):
        """RMSD between p_coords and q_coords. If weights are provided,
        RMSD is weighted."""
        diff = np.square(norm(p_coords - q_coords, axis=1))
        if weights is None:
            weights = 1
        diff = diff.reshape(-1,1)
        diff = weights*diff
        return np.sqrt(np.sum(diff) / diff.size)

    def _weights(self, p_coords, q_coords, c):
        """Calculate weights for each atom pair"""
        sq_diff = np.square(norm(p_coords - q_coords, axis=1))
        weights = np.exp((-sq_diff/c))
        weights = weights.reshape(-1,1)
        return weights

    def _auto_scaling_factor(self):
        """Automatic calculation of Gaussian scaling factor by getting the
        unweighted RMSD first, and feeding it in a sigmoid function"""
        self.run_unweighted()
        rms = self.rms
        scaling_factor = 8/(1+np.exp(-(1.7*rms-5.3)))+2
        return scaling_factor

    def _transform(self, coords, rot, tran):
        """Apply rotation matrix and translation vector to coordinates"""
        return dot(coords, rot) + tran

    def _fit(self, coords, reference_coords, weights=None):
        """Weighted superposition of two coordinate sets."""
        n = reference_coords.shape[0]
        if weights is None:
            weights = np.ones((n, 1))
        # Calculated weighted centroids
        centroid1 = np.sum(weights*coords, axis=0) / n
        centroid2 = np.sum(weights*reference_coords, axis=0) / n
        # Translate origin of both coordinates to centroid
        coords = coords - centroid1
        reference_coords = reference_coords - centroid2
        # Calculate weighted correlation matrix
        a = dot(transpose(coords)*weights.flatten(), reference_coords)
        # Singular value decomposition
        u, d, vt = svd(a)
        # Calculate rotation matrix
        rot = transpose(dot(transpose(vt), transpose(u)))
        # Check if matrix is a reflection
        if det(rot) < 0:
            vt[2] = -vt[2]
            rot = transpose(dot(transpose(vt), transpose(u)))
        # Calculate translation vector
        tran = centroid2 - dot(centroid1, rot)
        return rot, tran

