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

    for i, pdbsite in enumerate(entry.pdbsites):
        min = np.nanmean(pdbsite.get_distances('min'))
        ca = np.nanmean(pdbsite.get_distances('ca'))
        avg = np.nanmean(pdbsite.get_distances('avg'))
        com = np.nanmean(pdbsite.get_distances('com'))
        cog = np.nanmean(pdbsite.get_distances('cog'))
        print(pdbsite.id, min, ca, avg, com, cog)

        rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, weighted=True, scaling_factor=None, transform=True)
        ##per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite, rot, tran, transform=False)
        ##print(pdbsite.id, rms, rms_all, per_res_rms)
        pdbsite.write_pdb(outdir=outdir, write_hets=True, func_atoms_only=False, include_dummy_atoms=True)


if __name__ == '__main__':
    main(int(sys.argv[1]))
                
