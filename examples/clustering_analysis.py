#!/usr/bin/env python3
"""Cluster the catalytic sites of an M-CSA entry by 3D geometry.

Usage: python clustering_analysis.py [MCSA_ID]
"""

import sys

import numpy as np

import csa3d


def main(mcsa_id: int = 1) -> None:
    db = csa3d.Mcsa()  # assumes csa3d.setup_data() has been run
    db.build(entry_ids=mcsa_id)
    entry = db.entries[mcsa_id]

    sites = [site for site in entry.get_pdbsites(sane_only=True)
             if site.is_conserved or site.is_conservative_mutation]
    print(f"{len(sites)} conserved sites in M-CSA {mcsa_id}")

    matrix = entry.rmsd_matrix(sites)
    pairwise = matrix.values[np.triu_indices(len(matrix), k=1)]
    print(f"Pairwise RMSD: min {pairwise.min():.2f}, mean {pairwise.mean():.2f}, "
          f"max {pairwise.max():.2f} Å")

    print("\nHierarchical clustering (cut at 70% of tree height):")
    for cluster_id, site_ids in entry.clustering(matrix).items():
        print(f"  cluster {cluster_id}: {len(site_ids)} sites, e.g. {site_ids[0]}")

    print("\nBayesian cut of the same tree:")
    for cluster_id, site_ids in entry.clustering_bayesian(matrix).items():
        print(f"  cluster {cluster_id}: {len(site_ids)} sites")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
