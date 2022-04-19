#!/usr/bin/env python3

import sys
import os
import pickle
import Bio.PDB
import numpy as np
import pandas as pd

def main(mcsa_id):
    outdir = './out/mcsa_{}'.format(str(mcsa_id).zfill(4))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    print('Loading entry')
    with open('entries/csa3d_{}.ent'.format(str(mcsa_id).zfill(4)), 'rb') as f:
        entry = pickle.load(f)

    print('Building matrix')
    #matrix = entry.rmsd_matrix(entry.pdbsites)
    matrix = pd.read_csv('/Users/riziotis/ebi/phd/datasets/csa3d/per_entry_analyses/case_studies/csa3d_0133/rmsd_matrix.conserved.csv', index_col=0)

    print('Clustering')
    Z, clusters = entry.clustering(matrix)

    print('Building templates')
    for i, cluster in clusters.items():
        os.makedirs(f'{outdir}/cluster_{i}', exist_ok=True)
        print(f'Cluster {i}')
        for site in entry.pdbsites:
            if site.is_conserved or site.is_conservative_mutation and site.id in cluster:
                site.reference_site.fit(site, transform=True)
                site.write_pdb(outdir=f'{outdir}/cluster_{i}', func_atoms_only=True, write_hets=False)
        entry.create_template(outdir=outdir, subset=cluster, cluster_no=i)
    entry.create_template(outdir=outdir)

    
#    ref = entry.pdbsites[0].reference_site
#    for i, pdbsite in enumerate(entry.pdbsites):
#        rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, weighted=True, ca=False, scaling_factor=None, transform=True)
#    #    #per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite, rot, tran, transform=False)
#    #    ##print(pdbsite.id, rms, rms_all, per_res_rms)
#        pdbsite.write_pdb(outdir=outdir, write_hets=True, func_atoms_only=False, include_dummy_atoms=False)


if __name__ == '__main__':
    main(int(sys.argv[1]))
                
