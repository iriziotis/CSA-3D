#!/usr/bin/env python3

import pickle
import sys
sys.path.append('/Users/riziotis/ebi/phd/src')
sys.path.append('/nfs/research1/thornton/riziotis/research/phd/src')
from time import time
from mcsa3d import Mcsa

def main(mcsa_id, outdir):
    start = time()
    # Initialize database
    print('Initializing database')
    db = Mcsa()

    # Build test entry
    print('Building entry {}'.format(mcsa_id))
    i = time()
    build = db.build(mcsa_id, annotate=True, redundancy_cutoff=0.3, verbose=False)
    if not build:
        print('No entry {}, or has multiple references'.format(mcsa_id))
        exit()
    f = time()
    build_time = f-i

    # Serialiaze in a pickle obj and dump it in a file
    print('Serializing')
    i = time()
    with open('{}/csa3d_{}.ent'.format(outdir, str(mcsa_id).zfill(4)), 'wb') as o:
        pickle.dump(db.entries[mcsa_id], o)
    f = time()
    serialization_time = f-i

    finish = time()
    overall_time = finish-start

    print('Building time: {:5f}'.format(build_time))
    print('Serialization time: {:5f}'.format(serialization_time))
    print('Overall time: {:5f}'.format(overall_time))

if __name__ == '__main__':
    try:
        outdir = sys.argv[2]
    except IndexError:
        outdir = '.'
    main(int(sys.argv[1]), outdir)



