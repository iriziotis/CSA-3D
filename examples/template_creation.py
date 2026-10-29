#!/usr/bin/env python3
"""
CSA-3D Template Creation Example

This example demonstrates:
1. Loading and analyzing an entry
2. Creating consensus templates from site clusters
3. Breaking down templates into functional residue groups
"""


def template_creation():
    """Demonstrate template creation workflow."""
    try:
        from Mcsa import Mcsa
    except ImportError:
        print("Error: Could not import Mcsa module")
        print("Make sure you are running this from the CSA-3D root directory")
        return

    print("=" * 60)
    print("CSA-3D Template Creation Example")
    print("=" * 60)

    # Initialize database
    print("\n1. Initializing M-CSA database...")
    try:
        db = Mcsa()
        print("   ✓ Database initialized")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return

    # Build entry
    print("\n2. Building entry...")
    try:
        entry_id = 1
        success = db.build(entry_id, annotate=False)
        if not success:
            print(f"   ✗ Could not build entry {entry_id}")
            return

        entry = db.entries[entry_id]
        print(f"   ✓ Entry {entry_id} loaded")
        print(f"   - EC: {entry.e_number}")
        print(f"   - PDB sites: {len(entry.pdbsites)}")

    except Exception as e:
        print(f"   ✗ Error: {e}")
        return

    # Get conserved sites for clustering
    print("\n3. Preparing sites for clustering...")
    try:
        conserved = [s for s in entry.pdbsites if s.is_conserved]
        print(f"   ✓ Found {len(conserved)} conserved sites")

        if len(conserved) < 2:
            print("   - Need at least 2 sites for templates")
            return

    except Exception as e:
        print(f"   ✗ Error: {e}")
        return

    # Perform clustering
    print("\n4. Clustering sites...")
    try:
        matrix = entry.rmsd_matrix(conserved)
        clusters = entry.clustering_bayesian(matrix)
        print(f"   ✓ Created {len(clusters)} clusters")

    except Exception as e:
        print(f"   ✗ Error during clustering: {e}")
        return

    # Create templates for each cluster
    print("\n5. Creating templates...")
    try:
        for cluster_id, site_indices in clusters.items():
            print(f"\n   Cluster {cluster_id + 1}:")
            print(f"   - Sites in cluster: {len(site_indices)}")

            # Create template
            template = entry.create_template(
                ca=True,
                subset=site_indices,
                cluster_no=cluster_id + 1
            )
            print(f"   - Template created with {len(template)} residues")

            # Break template into functional groups
            try:
                groups = entry.break_template(
                    template,
                    n_residues=7,
                    permutations=False,
                    max_distance=8.0
                )
                print(f"   - Template divided into {len(groups)} residue groups")

                for i, (residues, atoms) in enumerate(groups.items(), 1):
                    print(f"     * Group {i}: {len(residues)} residues")

            except Exception as e:
                print(f"   - Template analysis: {e}")

    except Exception as e:
        print(f"   ✗ Error creating templates: {e}")
        return

    print("\n" + "=" * 60)
    print("Template creation completed!")
    print("\nTemplates can be used for:")
    print("  - Studying catalytic mechanism commonalities")
    print("  - Identifying key residues")
    print("  - Comparative analysis across clusters")
    print("=" * 60)


if __name__ == "__main__":
    template_creation()
