#!/usr/bin/env python3
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
import pickle

def main():
    with open('entry_2.ent', 'rb') as f:
        entry = pickle.load(f)

    print('Fitting')
    for pdbsite in entry.pdbsites:
        print(pdbsite.id, pdbsite)
        if pdbsite.is_reference:
            pdbsite.write_pdb(outdir='./out', write_hets=True)
            continue
        try:
            pdbsite.reference_site.fit(pdbsite, transform=True)
        except:
            continue
        pdbsite.write_pdb(outdir='./out', write_hets=True)



main()



