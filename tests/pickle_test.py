#!/usr/bin/env python3
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
import pickle

with open('entry_4.ent', 'rb') as f:
    entry = pickle.load(f)

entry.pdbsites[1].write_pdb(outdir='./test')
