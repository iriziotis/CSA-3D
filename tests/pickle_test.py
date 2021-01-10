#!/usr/bin/env python3
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
import pickle

def main():
    with open('entry_2.ent', 'rb') as f:
        entry = pickle.load(f)

    for pdbsite in entry.pdbsites:
        if len(pdbsite) > 7:
            continue
        if pdbsite.is_conserved or pdbsite.is_conservative_mutation:
            if len(pdbsite.mapped_unisites) == 0:
                print(pdbsite.pdb_id, pdbsite)
                continue
            for unisite in pdbsite.mapped_unisites:
                print(pdbsite.pdb_id, pdbsite, unisite.id)



main()


#print('Fitting')
#for pdbsite in entry.pdbsites:
#    if pdbsite.is_reference:
#        pdbsite.write_pdb(outdir='./out', write_hets=True)
#        continue
#    pdbsite.reference_site.fit(pdbsite, transform=True)
#    pdbsite.write_pdb(outdir='./out', write_hets=True)

