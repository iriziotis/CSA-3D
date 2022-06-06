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

    for s in entry.unisites:
        print(s.uniprot_id)
#    print('Building matrix')
#    sitelist = [site for site in entry.pdbsites if site.is_conserved or site.is_conservative_mutation]
#    matrix = entry.rmsd_matrix(sitelist)
#
#    print('Clustering')
#    clusters = entry.clustering_bayesian(matrix, plot_outfile=f'{outdir}/dendrogram.png')
#
#    print('Building templates')
#    for i, cluster in clusters.items():
#        os.makedirs(f'{outdir}/cluster_{i+1}', exist_ok=True)
#        print(f'Cluster {i}')
#
#        template = entry.create_template(ca=False, outdir=outdir, subset=cluster, cluster_no=i+1)
#
#       for site in entry.pdbsites:
#           if site.id in cluster:
#               site.reference_site.fit(site, transform=True)
#               site.write_pdb(outdir=f'{outdir}/cluster_{i+1}', func_atoms_only=True, write_hets=False)

    
#    ref = entry.pdbsites[0].reference_site
#    for i, pdbsite in enumerate(entry.pdbsites):
#        rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, weighted=True, ca=False, scaling_factor=None, transform=True)
#        per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite, rot, tran, transform=False)
#        print(pdbsite.id, rms, rms_all, per_res_rms)
#        pdbsite.write_pdb(outdir=outdir, write_hets=True, func_atoms_only=False, include_dummy_atoms=False)


if __name__ == '__main__':
    main(int(sys.argv[1]))
                
