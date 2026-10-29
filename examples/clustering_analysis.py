#!/usr/bin/env python3
"""
CSA-3D Clustering Analysis Example

This example demonstrates:
1. Loading an entry
2. Calculating RMSD between catalytic sites
3. Clustering sites by structure similarity
4. Visualizing clustering results
"""

import os
from pathlib import Path


def clustering_analysis():
    """Demonstrate RMSD clustering workflow."""
    try:
        from Mcsa import Mcsa
    except ImportError:
        print("Error: Could not import Mcsa module")
        print("Make sure you are running this from the CSA-3D root directory")
        return

    print("=" * 60)
    print("CSA-3D Clustering Analysis Example")
    print("=" * 60)

    # Initialize database
    print("\n1. Initializing M-CSA database...")
    try:
        db = Mcsa()
        print("   ✓ Database initialized")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return

    # Build entry with multiple sites
    print("\n2. Building entry with catalytic sites...")
    try:
        # Use entry with good coverage
        entry_id = 1
        success = db.build(entry_id, annotate=False)
        if not success:
            print(f"   ✗ Could not build entry {entry_id}")
            return

        entry = db.entries[entry_id]
        print(f"   ✓ Entry {entry_id} loaded")
        print(f"   - EC number: {entry.e_number}")
        print(f"   - Total PDB sites: {len(entry.pdbsites)}")

    except Exception as e:
        print(f"   ✗ Error building entry: {e}")
        return

    # Filter sites (conserved or conservative mutations)
    print("\n3. Filtering conserved sites...")
    try:
        conserved_sites = [
            site for site in entry.pdbsites
            if site.is_conserved or site.is_conservative_mutation
        ]
        print(f"   ✓ Found {len(conserved_sites)} conserved/conservative sites")

        if len(conserved_sites) < 2:
            print("   - Not enough sites for clustering (need >= 2)")
            return
    except Exception as e:
        print(f"   ✗ Error filtering sites: {e}")
        return

    # Calculate RMSD matrix
    print("\n4. Calculating RMSD matrix...")
    try:
        rmsd_matrix = entry.rmsd_matrix(conserved_sites)
        print(f"   ✓ RMSD matrix calculated")
        print(f"   - Shape: {rmsd_matrix.shape}")
        print(f"   - Min: {rmsd_matrix.min():.3f} Å")
        print(f"   - Max: {rmsd_matrix.max():.3f} Å")
        print(f"   - Mean: {rmsd_matrix.mean():.3f} Å")

    except Exception as e:
        print(f"   ✗ Error calculating RMSD: {e}")
        return

    # Perform clustering
    print("\n5. Performing hierarchical clustering...")
    try:
        clusters = entry.clustering_bayesian(rmsd_matrix)
        print(f"   ✓ Clustering completed")
        print(f"   - Number of clusters: {len(clusters)}")

        for cluster_id, site_indices in clusters.items():
            print(f"   - Cluster {cluster_id + 1}: {len(site_indices)} sites")

    except Exception as e:
        print(f"   ✗ Error during clustering: {e}")
        print("   Note: Clustering requires external dependencies")
        return

    # Summary
    print("\n" + "=" * 60)
    print("Clustering analysis completed!")
    print("\nTip: Use the clustering results to:")
    print("  - Identify structurally similar catalytic sites")
    print("  - Create templates for each cluster")
    print("  - Study catalytic mechanism variations")
    print("=" * 60)


if __name__ == "__main__":
    clustering_analysis()
