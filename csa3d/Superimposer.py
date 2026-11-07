"""Weighted structural superposition for catalytic site alignment.

Implements hybrid iterative and weighted superposition to align protein
structures accounting for dynamic conformational variation. The algorithm
uses outlier rejection followed by Gaussian-weighted fitting to handle
geometric diversity in catalytic sites.

Reference:
    Riziotis et al., Journal of Molecular Biology, 2022.
    "Conformational Variation in Enzyme Catalysis: A Structural Study on 
    Catalytic Residues" https://doi.org/10.1016/j.jmb.2022.167517
    
    Also see: Coutsias et al., Journal of Computational Chemistry, 2004.
    "Using quaternions to calculate RMSD"
    https://dx.doi.org/10.1529%2Fbiophysj.105.066654
"""

import numpy as np
from numpy import dot, transpose, sqrt
from numpy.linalg import svd, det, norm
from copy import copy


class Superimposer:
    """Superimpose two sets of coordinates using weighted least-squares fitting.
    
    Performs structural superposition of coordinate sets using either:
    1. Iterative unweighted fitting with outlier rejection (fast, robust to outliers)
    2. Hybrid weighted fitting combining both approaches (precise)
    
    The unweighted method iteratively superimposes coordinates, removing atoms
    beyond a distance cutoff at each iteration. This is fast and handles large
    conformational deviations. The weighted method follows this with Gaussian-weighted
    iterations to achieve sub-Ångstrom precision by downweighting atoms in
    regions of uncertainty.
    
    This approach is particularly useful for aligning catalytic sites where residue
    conformations can vary significantly across homologous structures, while the
    overall active site geometry remains conserved. The Gaussian weighting function
    emphasizes well-aligned atoms and de-emphasizes atoms in flexible regions.
    
    Attributes:
        reference_coords: Reference coordinate set (N×3 NumPy array)
        coords: Mobile coordinate set to be superimposed (N×3 array)
        rot: Optimal rotation matrix (3×3)
        tran: Optimal translation vector (3,)
        rms: Weighted RMSD after fitting (excludes outliers if used)
        rms_all: Unweighted RMSD over all atoms
        cycles: Number of outlier rejection iterations
        cutoff: Distance threshold for outlier rejection (Ångströms)
        scaling_factor: Gaussian weighting function scaling parameter
        n: Number of atoms
    """

    def __init__(self):
        """Initialize a Superimposer instance with cleared state."""
        self._clear()

    # Public methods

    def set(self, reference_coords, coords, cycles=1, cutoff=999, scaling_factor=None):
        """Set coordinates and parameters for superposition.
        
        Configures the coordinate sets and fitting parameters. Automatically computes
        an adaptive Gaussian scaling factor if not provided. Must be called before
        run_unweighted() or run_weighted().
        
        Args:
            reference_coords: Reference coordinate set as N×3 NumPy array,
                where N is the number of points and columns are (x, y, z).
            coords: Mobile coordinate set to be superimposed, must have same shape.
            cycles: Number of iterations of iterative unweighted superposition 
                with outlier rejection. Higher values are more robust to outliers 
                but slower. Defaults to 1. Set to >1 for structures with large 
                conformational deviations.
            cutoff: Distance threshold (Ångströms) for outlier rejection during
                iterative fitting. Atoms with inter-coordinate distance >cutoff
                are excluded in the next iteration. Defaults to 999 (essentially 
                no cutoff). Typical values are 2-6 Å depending on expected flexibility.
            scaling_factor: Gaussian scaling parameter (σ²) for weighted fitting.
                Controls how much weight is given to well-aligned atoms. If None
                (default), computed automatically from the unweighted RMSD using
                a sigmoid function. Manual specification useful for reproducible
                results or fine-tuning behavior.
        
        Raises:
            Exception: If reference_coords and coords have different dimensions.
        
        Example:
            >>> import numpy as np
            >>> sup = Superimposer()
            >>> ref = np.random.rand(10, 3)
            >>> mob = np.random.rand(10, 3)
            >>> sup.set(ref, mob, cycles=3, cutoff=3.0)
        """
        self._clear()
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
        """Iteratively superimpose coordinates with outlier rejection.
        
        Performs classic rigid-body superposition using the Kabsch algorithm
        (SVD-based least-squares fitting). At each iteration, atoms with
        inter-coordinate distances exceeding the cutoff are removed, and fitting
        is repeated on the remaining atoms. This is fast and robust to outliers
        but may not achieve sub-Ångstrom precision if significant conformational
        variation is present.
        
        Results are stored in self.rot, self.tran, self.rms, self.rms_all.
        Call get_rotran() or get_rms() to retrieve them.
        
        The method maintains at least 3 atoms for superposition (minimum for
        rigid-body transformation). If fewer than 3 atoms remain after outlier
        rejection, iteration stops early.
        
        Note:
            Must call set() first. Results are automatically used by run_weighted()
            as the initialization step.
        """
        reference_coords = self.reference_coords.copy()
        coords = self.coords.copy()
        coords_all = self.coords.copy()
        untransformed_coords = self.coords.copy()
        min_rms = 999

        # Run superposition
        for i in range(self.cycles):
            # Fit using SVD-based least squares
            rot, tran = self._fit(coords, reference_coords)
            # Transform coordinates
            coords = self._transform(coords, rot, tran)
            coords_all = self._transform(coords_all, rot, tran)
            # Calculate RMSDs
            rms = self._rms(reference_coords, coords)
            rms_all = self._rms(self.reference_coords, coords_all)
            # Reject outliers beyond cutoff
            if rms < min_rms:
                min_rms = rms
            diff = norm(reference_coords-coords, axis=1)
            to_keep = np.where(diff < self.cutoff)
            # Maintain minimum of 3 atoms for superposition
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
        """Hybrid iterative and Gaussian-weighted superposition.
        
        Performs a two-stage fitting process:
        1. Unweighted iterative fitting with outlier rejection (calls run_unweighted())
        2. Gaussian-weighted iterative fitting to convergence
        
        The second stage uses a Gaussian weighting function where weight decreases
        with squared distance. This de-emphasizes atoms in regions of conformational
        uncertainty while emphasizing the conserved core geometry. Iteration continues
        until the RMSD improvement between cycles is <0.0001 Ångströms.
        
        Results are stored in self.rot, self.tran, self.rms, self.rms_all.
        Call get_rotran() or get_rms() to retrieve them.
        
        The scaling_factor (σ² in Gaussian weight) controls how sharply atoms
        are down-weighted by distance. Larger values give more weight to distant
        atoms (closer to unweighted fitting). Smaller values emphasize only the
        best-aligned atoms. The default auto-computed value is often appropriate,
        but can be overridden in set() for specific applications.
        
        Raises:
            Exception: If convergence is not achieved after 6000 weighted iterations.
        
        Note:
            Must call set() first. This method is slower but more accurate than
            run_unweighted() for catalytic sites with multiple conformers.
        """
        coords_all = self.coords.copy()
        n = self.n
        c = self.scaling_factor

        # Initial iterative superposition
        self.run_unweighted()
        coords = self.temp_coords.copy()
        reference_coords = self.temp_reference_coords.copy()
        untransformed_coords = self.temp_untransformed_coords.copy()
        coords_all = self._transform(coords_all, self.rot, self.tran)

        # Weighted iterative superposition until convergence
        prev_rms = 999
        for i in range(6000):
            # Compute Gaussian weights based on current distances
            weights = self._weights(coords, reference_coords, c)
            # Perform weighted fit
            rot, tran = self._fit(coords, reference_coords, weights)
            # Transform coordinates
            coords = self._transform(coords, rot, tran)
            coords_all = self._transform(coords_all, rot, tran)
            # Calculate weighted RMSD 
            rms = self._rms(coords, reference_coords, weights)
            # RMSD over all atoms, unweighted
            rms_all = self._rms(self.reference_coords, coords_all)
            # Check for convergence
            if 0 <= prev_rms - rms <= 0.0001:
                break
            prev_rms = rms

        if i == 5999:
            raise Exception("Could not reach convergence after 6000 cycles. Maybe try unweighted fit.")
        else:
            # Final rotation matrix, translation vector and RMSD
            self.rot, self.tran = self._fit(untransformed_coords, coords)
            self.rms = rms
            self.rms_all = rms_all

    def get_rotran(self):
        """Return the rotation matrix and translation vector.
        
        These represent the transformation applied to the mobile coordinates
        to best superimpose them onto the reference. The transformation is:
        
            transformed_coords = mobile_coords @ rot + tran
        
        Returns:
            tuple: (rot, tran) where:
                - rot (3×3 ndarray): Rotation matrix (right multiplying)
                - tran (3,) ndarray: Translation vector
        
        Raises:
            Exception: If superposition has not been performed yet.
        
        Example:
            >>> sup = Superimposer()
            >>> sup.set(ref_coords, mob_coords)
            >>> sup.run_weighted()
            >>> rot, tran = sup.get_rotran()
            >>> transformed = mob_coords @ rot + tran
        """
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rot, self.tran

    def get_rms(self):
        """Return the RMSD values from the last superposition.
        
        Returns two RMSD measures:
            - rms: RMSD after outlier rejection (if used), excluding atoms rejected
          as outliers during iterative fitting. For unweighted fitting with
          outlier rejection, this is the RMSD of the final subset of atoms.
            - rms_all: Unweighted RMSD over all original atoms, regardless of whether
          they were used in the iterative fitting.
        
        For weighted fitting, rms is the weighted RMSD and rms_all is computed
        over all atoms with equal weight.
        
        Returns:
            tuple: (rms, rms_all) where both are floats (Ångströms)
        
        Raises:
            Exception: If superposition has not been performed yet.
        
        Example:
            >>> sup = Superimposer()
            >>> sup.set(ref_coords, mob_coords)
            >>> sup.run_unweighted()
            >>> rms, rms_all = sup.get_rms()
            >>> print(f"Core RMSD: {rms:.2f} Å, Full RMSD: {rms_all:.2f} Å")
        """
        if self.rot is None:
            raise Exception("Nothing superimposed yet.")
        return self.rms, self.rms_all

    # Private methods

    def _clear(self):
        """Reset all internal state to None."""
        self.reference_coords = None
        self.coords = None
        self.rot = None
        self.tran = None
        self.rms = None
        self.cycles = None
        self.cutoff = None
        self.scaling_factor = None

    def _rms(self, p_coords, q_coords, weights=None):
        """Compute RMSD between two coordinate sets.
        
        If weights are provided, computes the weighted RMSD where each atom
        contributes to the sum proportionally to its weight. Otherwise, all
        atoms contribute equally.
        
        Args:
            p_coords: First coordinate set (N×3 array)
            q_coords: Second coordinate set (N×3 array)
            weights: Optional weights (N×1 array). Defaults to None (equal weights).
        
        Returns:
            float: RMSD in Ångströms
        """
        diff = np.square(norm(p_coords - q_coords, axis=1))
        if weights is None:
            weights = 1
            n = diff.size
        else:
            n = np.sum(weights)
        diff = diff.reshape(-1, 1)
        diff = weights * diff
        return np.sqrt(np.sum(diff) / n)

    def _weights(self, p_coords, q_coords, c):
        """Compute Gaussian weights based on distance between atoms.
        
        Weight = exp(-(squared_distance / c))
        
        This gives weight ≈ 1 to well-aligned atoms and exponentially lower
        weights to atoms with larger distances. The parameter c controls the
        decay rate: smaller c emphasizes good alignment more sharply.
        
        Args:
            p_coords: Mobile coordinates after transformation (N×3 array)
            q_coords: Reference coordinates (N×3 array)
            c: Gaussian scaling parameter (σ²). Typically 1-20 depending on
               expected flexibility and RMSD.
        
        Returns:
            ndarray: Weights (N×1), each in range [0, 1]
        """
        sq_diff = np.square(norm(p_coords - q_coords, axis=1))
        weights = np.exp((-sq_diff / c))
        weights = weights.reshape(-1, 1)
        return weights

    def _auto_scaling_factor(self):
        """Compute adaptive Gaussian scaling parameter from data.
        
        Runs unweighted superposition and uses its RMSD to compute an
        adaptive scaling factor via a sigmoid function. This makes the
        weighting sensitive to the initial RMSD: smaller initial RMSD
        gives smaller scaling factor (sharper weighting), while large
        RMSD gives larger scaling factor (softer weighting).
        
        The sigmoid formula is: σ² = 9/(1 + exp(-(1.7·RMSD - 5.3))) + 1
        This produces scaling factors typically in the range [1, 10].
        
        Returns:
            float: Adaptive scaling factor for Gaussian weighting
        """
        self.run_unweighted()
        rms = self.rms
        scaling_factor = 9 / (1 + np.exp(-(1.7 * rms - 5.3))) + 1
        return scaling_factor

    def _transform(self, coords, rot, tran):
        """Apply rotation and translation to coordinates.
        
        Args:
            coords: Coordinate array (N×3)
            rot: Rotation matrix (3×3)
            tran: Translation vector (3,)
        
        Returns:
            ndarray: Transformed coordinates (N×3)
        """
        return dot(coords, rot) + tran

    def _fit(self, coords, reference_coords, weights=None):
        """Compute optimal rotation and translation using weighted least squares.
        
        Implements the Kabsch algorithm with optional weighting. Finds the
        rotation matrix and translation vector that minimizes the weighted
        sum of squared distances.
        
        Algorithm:
        1. Compute weighted centroids of both sets
        2. Center both around their centroids
        3. Compute weighted correlation matrix (covariance-like)
        4. SVD of correlation matrix
        5. Extract rotation matrix (ensuring proper orientation, not reflection)
        6. Compute translation from centroids
        
        Args:
            coords: Mobile coordinates (N×3 array)
            reference_coords: Reference coordinates (N×3 array)
            weights: Optional per-atom weights (N×1 array). Defaults to None (equal).
        
        Returns:
            tuple: (rot, tran) - rotation matrix (3×3) and translation vector (3,)
        """
        n = reference_coords.shape[0]
        if weights is None:
            weights = np.ones((n, 1))
        # Compute weighted centroids
        centroid1 = np.sum(weights * coords, axis=0) / n
        centroid2 = np.sum(weights * reference_coords, axis=0) / n
        # Translate to origin
        coords = coords - centroid1
        reference_coords = reference_coords - centroid2
        # Compute weighted correlation matrix
        a = dot(transpose(coords) * weights.flatten(), reference_coords)
        # SVD decomposition
        u, d, vt = svd(a)
        # Construct rotation matrix
        rot = transpose(dot(transpose(vt), transpose(u)))
        # Ensure proper rotation (not reflection)
        if det(rot) < 0:
            vt[2] = -vt[2]
            rot = transpose(dot(transpose(vt), transpose(u)))
        # Compute translation vector
        tran = centroid2 - dot(centroid1, rot)
        return rot, tran
