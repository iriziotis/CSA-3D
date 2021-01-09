#!/usr/bin/env python3
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
import pickle

print('Loading entry')
with open('entry_1.ent', 'rb') as f:
    entry = pickle.load(f)

print('Fitting')
for pdbsite in entry.pdbsites:
    if pdbsite.is_reference:
        pdbsite.write_pdb(outdir='./out', write_hets=True)
        continue
    pdbsite.reference_site.fit(pdbsite, transform=True)
    pdbsite.write_pdb(outdir='./out', write_hets=True)

