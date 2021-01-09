#!/usr/bin/env python3
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
import pickle

with open('entry_1.ent', 'rb') as f:
    entry = pickle.load(f)

#entry.pdbsites[1].write_pdb(outdir='./out')
for pdbsite in entry.pdbsites:
    for i in pdbsite.mapped_unisites:
        print(pdbsite.id, i.id)
