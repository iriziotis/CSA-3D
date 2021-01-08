#!/usr/bin/env python3

import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
from mcsa3d import Mcsa
import pickle

print('Initializing database')
db = Mcsa()
print('Building')
db.build([1], annotate=False, verbose=True)

#print('Fitting')
#for site in db.entries[1].pdbsites:
#    if site.is_reference:
#        site.write_pdb(outdir='./test', write_hets=True)
#        continue
#    try:
#        site.reference_site.fit(site, transform=True, reorder=False, mutate=False, cycles=5)
#    except Exception as e:
#        print(site.id, e)
#        continue
#    site.write_pdb(outdir='./test', write_hets=True)

print('Serializing')
with open('entry_4.ent', 'wb') as o:
    pickle.dump(db.entries[1], o)



    
    

