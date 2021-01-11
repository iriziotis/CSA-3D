#!/usr/bin/env python3
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
import os
import pickle

def main(mcsa_id):
    with open('entry_{}.ent'.format(mcsa_id), 'rb') as f:
        entry = pickle.load(f)

    outdir = './out/mcsa_{}'.format(str(mcsa_id).zfill(4))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for pdbsite in entry.pdbsites:
        rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, transform=True)
        print(pdbsite.id, rms, rms_all)
        pdbsite.write_pdb(outdir=outdir, write_hets=True, func_atoms_only=False)

if __name__ == '__main__':
    main(int(sys.argv[1]))

