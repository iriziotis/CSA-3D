#!/usr/bin/env python3

import sys
import os
import pickle
import Bio.PDB

def main(mcsa_id):
    with open('entries/csa3d_{}.ent'.format(str(mcsa_id).zfill(4)), 'rb') as f:
        entry = pickle.load(f)

    outdir = './out/mcsa_{}'.format(str(mcsa_id).zfill(4))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    seen = set()
    for i in entry.pdbsites:
        for j in entry.pdbsites:
            if i==j or (j.id, i.id) in seen:
                continue
            seen.add((i.id, j.id))
            r, t, rms, rms_all = i.fit(j, weighted=True)
            pr = i.per_residue_rms(j)
            print(i.id, j.id, rms, rms_all, pr, max(pr), max(pr)<rms_all)

    #for pdbsite in entry.pdbsites:
    #    rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, weighted=True, cycles=10, cutoff=6, scaling_factor=None, transform=True)
    #    per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite)
    #    print(pdbsite.id, rms, rms_all)
    #    pdbsite.write_pdb(outdir=outdir, write_hets=False, func_atoms_only=True)

if __name__ == '__main__':
    main(int(sys.argv[1]))
                
