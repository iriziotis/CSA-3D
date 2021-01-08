#!/usr/bin/env python3

import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
from mcsa3d import Mcsa

print('Initializing database')
db = Mcsa()
print('Building')
db.build([4], annotate=False, verbose=True)

print('Fitting')
for site in db.entries[4].pdbsites:
    if site.is_reference:
        site.write_pdb(outdir='./test', write_hets=True)
        continue
    try:
        site.reference_site.fit(site, transform=True, reorder=False, mutate=False, cycles=1)
    except Exception as e:
        print(site.id, e)
    site.write_pdb(outdir='./test', write_hets=True)



    
    

