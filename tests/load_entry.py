#!/usr/bin/env python3

import sys
import os
import pickle
import Bio.PDB
import numpy as np

def main(mcsa_id):
    with open('entries/csa3d_{}.ent'.format(str(mcsa_id).zfill(4)), 'rb') as f:
        entry = pickle.load(f)

    outdir = './out/mcsa_{}'.format(str(mcsa_id).zfill(4))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    sitelist = [s for s in entry.get_pdbsites(sane_only=True) if  (s.is_conserved or s.is_conservative_mutation)]
    matrix = entry.rmsd_matrix(sitelist)
    Z, clusters = entry.clustering(matrix)
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    fig, ax = plt.subplots()
    dendrogram(Z, labels=matrix.index, ax=ax)

    print(clusters)
    plt.show()



    #for i, pdbsite in enumerate(entry.pdbsites):
    #    rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, weighted=True, scaling_factor=None, transform=True)
    #    #per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite, rot, tran, transform=False)
    #    #print(pdbsite.id, rms, rms_all, per_res_rms)
    #    pdbsite.write_pdb(outdir=outdir, write_hets=True, func_atoms_only=False, include_dummy_atoms=True)


if __name__ == '__main__':
    main(int(sys.argv[1]))
                
