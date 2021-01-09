#!/usr/bin/env python3

import pickle
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
from mcsa3d import Mcsa

entry = 1

print('Initializing database')
db = Mcsa()
print('Building')
db.build([entry], annotate=False, verbose=True)

print('Fitting')
for site in db.entries[entry].pdbsites:
    if site.is_reference:
        site.write_pdb(outdir='./out', write_hets=True)
        continue
    try:
        site.reference_site.fit(site, transform=True, reorder=False, mutate=False, cycles=5)
    except Exception as e:
        continue
    site.write_pdb(outdir='./out', write_hets=True)

#print('Serializing')
#with open('entry_{}.ent'.format(entry), 'wb') as o:
#    pickle.dump(db.entries[entry], o)



    
    

