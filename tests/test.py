#!/usr/bin/env python3

import pickle
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
from time import time
from mcsa3d import Mcsa

entry = 4

# Initialize database
print('Initializing database')
db = Mcsa()

# Build test entry
print('Building')
start = time()
db.build([entry], annotate=False, verbose=False)
finish = time()
build_time = finish-start

# Serialiaze in a pickle obj and dump it in a file
print('Serializing')
start = time()
with open('entry_{}.ent'.format(entry), 'wb') as o:
    pickle.dump(db.entries[entry], o)
finish = time()
serialization_time = finish-start

print('Building time: {:3f}'.format(build_time))
print('Serialization time: {:3f}'.format(serialization_time))

#print('Fitting')
#for site in db.entries[entry].pdbsites:
#    if site.is_reference:
#        site.write_pdb(outdir='./out', write_hets=True)
#        continue
#    try:
#        site.reference_site.fit(site, transform=True, reorder=False, mutate=False, cycles=5)
#    except Exception as e:
#        continue
#    site.write_pdb(outdir='./out', write_hets=True)




    
    

