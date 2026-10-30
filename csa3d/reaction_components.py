#!/usr/bin/env python3
"""Makes a mapping csv file for rhea reaction ids to chebi components"""

import csv
from collections import defaultdict
from natsort import natsorted
from glob import glob

def main():

    # Map directional reaction ids to master ids
    DIRECTIONAL2MASTER = {}
    with open('rhea-directions.tsv', 'r') as f:
        next(f)
        for line in csv.reader(f, delimiter='\t'):
            master = line[0]
            lr = line[1]
            rl = line[2]
            DIRECTIONAL2MASTER[lr] = master
            DIRECTIONAL2MASTER[rl] = master

    # Map ECs to master IDs
    MASTER2EC = defaultdict(set)
    with open('./rhea2ec.tsv', 'r') as f:
        next(f)
        for line in csv.reader(f, delimiter='\t'):
            master = line[2]
            ec = line[3]
            MASTER2EC[master].add(ec)

    # Make mapping csv
    REACTIONS = defaultdict(set)
    for infile in natsorted(glob('rd/*.rd')):
        rhea_id = infile.split('/')[-1].split('.')[0]
        with open(infile, 'r') as f:
            for line in f:
                if line.startswith('CHEBI:'):
                    chebi_id = line.strip().split(':')[-1]
                    # Test if something is wrong
                    int(rhea_id), int(chebi_id)
                    REACTIONS[rhea_id].add(chebi_id)

    # Print to stdout
    print('RHEA_MASTER_ID,RHEA_DIRECTIONAL_ID,EC,CHEBI_IDS')
    for k,v in REACTIONS.items():
        master = DIRECTIONAL2MASTER.get(k, int(k)-1)
        ecs = MASTER2EC.get(master,set())
        for ec in ecs:
            print('{},{},{},{}'.format(master, k, ec, ';'.join(v)))

if __name__ == '__main__':
    main()
            
