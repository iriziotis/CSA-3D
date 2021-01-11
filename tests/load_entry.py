#!/usr/bin/env python3
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
import pickle
import numpy as np

def main(mcsa_id):
    with open('entry_{}.ent'.format(mcsa_id), 'rb') as f:
        entry = pickle.load(f)

    for pdbsite in entry.pdbsites:
        print(pdbsite.id, pdbsite.uniprot_id, ','.join(i.id for i in pdbsite.mapped_unisites))
        #rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, transform=True)
        #print(pdbsite.id, rms, rms_all)
        #pdbsite.write_pdb(outdir='./out', write_hets=True)

if __name__ == '__main__':
    main(int(sys.argv[1]))

