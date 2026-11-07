#!/usr/bin/env python3
"""Basic CSA-3D workflow: build an M-CSA entry and inspect its catalytic sites.

Usage: python basic_usage.py [MCSA_ID]
"""

import sys

import csa3d


def main(mcsa_id: int = 1) -> None:
    csa3d.setup_data()  # one-off download; skips files already present

    db = csa3d.Mcsa()
    db.build(entry_ids=mcsa_id)
    entry = db.entries[mcsa_id]
    print(f"M-CSA {mcsa_id}: {entry.info['enzyme_name']} (EC {', '.join(entry.info['all_ecs'])})")
    print(f"{len(entry.pdbsites)} PDB sites, {len(entry.unisites)} UniProt sites\n")

    reference = next(site for site in entry.get_pdbsites() if site.is_reference)
    print(f"Reference site {reference.id} ({reference.resolution} Å)")
    for res in reference.residues:
        print(f"  {res.resname} {res.chain}{res.auth_resid}")

    print("\nHomologous sites (first 5):")
    sites = [site for site in entry.get_pdbsites(sane_only=True) if not site.is_reference]
    for site in sites[:5]:
        _, _, rms, _ = reference.fit(site)
        status = "conserved" if site.is_conserved else "mutated"
        print(f"  {site.id:40s} {site.sequence}  {status:9s} RMSD {rms:.2f} Å")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
