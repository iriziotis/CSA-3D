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

    ## Entry 789
    #a = entry.get_pdbsite('5j7t_A-A-A-A')
    #b = entry.get_pdbsite('5chv_A-A-A-A')
    #_,_,_, rb = b.fit(a, weighted=True, scaling_factor=None, transform=True)
    #print(rb)
    #a.write_pdb(outfile='a.pdb', func_atoms_only=True)
    #b.write_pdb(outfile='b.pdb', func_atoms_only=True)

    for i, pdbsite in enumerate(entry.pdbsites):
        rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, weighted=True, cycles=10, cutoff=6, scaling_factor=None, transform=True)
        per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite)
        print(pdbsite.id, rms, rms_all, per_res_rms)
        pdbsite.write_pdb(outdir=outdir, write_hets=False, func_atoms_only=True)


if __name__ == '__main__':
    main(int(sys.argv[1]))
                
