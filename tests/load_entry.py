#!/usr/bin/env python3

import sys
import os
import pickle

def main(mcsa_id):
    with open('csa3d_{}.ent'.format(str(mcsa_id).zfill(4)), 'rb') as f:
        entry = pickle.load(f)

    outdir = './out/mcsa_{}'.format(str(mcsa_id).zfill(4))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for pdbsite in entry.pdbsites:
        rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, transform=True)
        per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite)
        print(pdbsite.id, per_res_rms)
        #print(pdbsite.id, pdbsite.sequence, rms, rms_all, ','.join(str(i) for i in pdbsite.mapped_unisites))
        pdbsite.write_pdb(outdir=outdir, write_hets=False, func_atoms_only=True)
    
if __name__ == '__main__':
    main(int(sys.argv[1]))

