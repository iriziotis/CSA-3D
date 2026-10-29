#!/usr/bin/env python3
"""
Basic CSA-3D Usage Example

This example demonstrates the fundamental workflow for using CSA-3D:
1. Initializing the database
2. Building entries
3. Accessing catalytic sites
"""

# Note: This example requires the M-CSA data files to be available.
# Set the CSA3D_DATA_DIR environment variable to point to your data directory,
# or place data in ~/.csa3d/data

def basic_workflow():
    """Demonstrate basic CSA-3D workflow."""
    try:
        from Mcsa import Mcsa
    except ImportError:
        print("Error: Could not import Mcsa module")
        print("Make sure you are running this from the CSA-3D root directory")
        return

    print("=" * 60)
    print("CSA-3D Basic Usage Example")
    print("=" * 60)

    # Initialize the database
    print("\n1. Initializing M-CSA database...")
    try:
        db = Mcsa()
        print(f"   ✓ Database initialized")
    except Exception as e:
        print(f"   ✗ Error initializing database: {e}")
        print(f"   Make sure CSA3D_DATA_DIR is set correctly")
        return

    # Build a specific entry
    print("\n2. Building an entry (ID: 1)...")
    try:
        success = db.build(1, annotate=False)
        if success:
            entry = db.entries[1]
            print(f"   ✓ Entry 1 built successfully")
            print(f"   - EC number: {entry.e_number}")
            print(f"   - PDB sites: {len(entry.pdbsites)}")
            print(f"   - UniProt sites: {len(entry.unisites)}")
        else:
            print(f"   ✗ Failed to build entry 1")
            return
    except Exception as e:
        print(f"   ✗ Error building entry: {e}")
        return

    # Access catalytic sites
    print("\n3. Accessing catalytic sites...")
    try:
        for i, site in enumerate(entry.pdbsites[:3], 1):
            print(f"   PDB Site {i}:")
            print(f"   - ID: {site.id}")
            print(f"   - Structure: {site.pdb_id}")
            print(f"   - Residues: {len(site)}")
            if i >= 3:
                break

        if len(entry.unisites) > 0:
            print(f"\n   UniProt Sites: {len(entry.unisites)} total")
    except Exception as e:
        print(f"   ✗ Error accessing sites: {e}")
        return

    # Calculate RMSD matrix (if multiple sites)
    print("\n4. Calculating RMSD matrix...")
    try:
        if len(entry.pdbsites) > 1:
            sites = entry.pdbsites[:5]  # Use first 5 sites
            matrix = entry.rmsd_matrix(sites)
            print(f"   ✓ RMSD matrix calculated")
            print(f"   - Matrix shape: {matrix.shape}")
            print(f"   - Matrix range: {matrix.min():.2f} to {matrix.max():.2f}")
        else:
            print(f"   - Only {len(entry.pdbsites)} site(s) available")
    except Exception as e:
        print(f"   ✗ Error calculating RMSD: {e}")

    print("\n" + "=" * 60)
    print("Example completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    basic_workflow()
