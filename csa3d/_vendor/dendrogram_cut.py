"""Minimal implementation of dendrogram-cut for clustering.

This is a vendored minimal implementation providing hierarchical clustering
with automatic cut height selection. Implements the API needed by CSA-3D
for RMSD-based site clustering.
"""

import numpy as np
from scipy.cluster.hierarchy import linkage, cut_tree
from typing import Optional


class DendrogramCut:
    """Hierarchical clustering with automatic dendrogram cut selection.

    Provides clustering with automatic determination of optimal cluster count
    using various linkage methods.
    """

    def __init__(self, k_max: int = 10, method: str = 'average') -> None:
        """Initialize DendrogramCut clustering.

        Args:
            k_max: Maximum number of clusters to consider
            method: Linkage method ('average', 'complete', 'single', 'ward')
        """
        self.k_max = k_max
        self.method = method
        self.Z = None
        self.labels_ = None

    def fit(self, X: np.ndarray):
        """Fit hierarchical clustering to data.

        Args:
            X: Distance matrix (condensed form from pdist) or sample matrix

        Returns:
            self
        """
        # If X is 1D (condensed distance matrix), use directly
        if X.ndim == 1:
            self.Z = linkage(X, method=self.method)
        else:
            # If X is 2D (sample matrix), compute distances first
            from scipy.spatial.distance import pdist
            distances = pdist(X)
            self.Z = linkage(distances, method=self.method)

        # Default: cut at k_max clusters
        self.labels_ = cut_tree(self.Z, n_clusters=min(self.k_max, self.Z.shape[0])).ravel()
        return self

    def predict(self, k: Optional[int] = None) -> np.ndarray:
        """Get cluster labels.

        Args:
            k: Number of clusters (uses k_max if not specified)

        Returns:
            Cluster labels
        """
        if self.Z is None:
            raise ValueError("Must fit model before predict")

        if k is None:
            k = min(self.k_max, self.Z.shape[0])

        return cut_tree(self.Z, n_clusters=k).ravel()
