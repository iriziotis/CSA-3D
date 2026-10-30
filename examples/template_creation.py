#!/usr/bin/env python3
"""Build 3D templates of an M-CSA catalytic site for Jess / PyJess.

One template is written per geometric cluster of sites, plus smaller templates of
3 neighbouring residues each. Templates can be searched against protein structures
with Jess (https://github.com/iriziotis/jess) or PyJess (https://github.com/althonos/pyjess).

Usage: python template_creation.py [MCSA_ID] [OUTDIR]
"""

import sys
from pathlib import Path

import csa3d


def main(mcsa_id: int = 1, outdir: str = "templates") -> None:
    Path(outdir).mkdir(exist_ok=True)
    db = csa3d.Mcsa()  # assumes csa3d.setup_data() has been run
    db.build(entry_ids=mcsa_id)
    entry = db.entries[mcsa_id]

    sites = list(entry.get_pdbsites(sane_only=True))
    clusters = entry.clustering(entry.rmsd_matrix(sites))

    for cluster_id, site_ids in clusters.items():
        template = entry.create_template(outdir=outdir, subset=site_ids, cluster_no=cluster_id)
        if template is None:
            print(f"Cluster {cluster_id}: no conserved sites to build a template from")
            continue
        print(f"Cluster {cluster_id}: template of {template.size} residues from {len(site_ids)} sites")

        for group in entry.break_template(template, n_residues=3):
            entry.write_template(template, residues=list(group), subset=site_ids,
                                 cluster_no=cluster_id, outdir=outdir)

    print(f"\nWritten to {outdir}/:")
    for path in sorted(Path(outdir).iterdir()):
        print(" ", path.name)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 1, *sys.argv[2:3])
