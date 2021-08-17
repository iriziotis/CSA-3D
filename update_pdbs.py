#!/usr/bin/env python3

from csa3d import Mcsa

print('Loading data...')
db = Mcsa()

print('Downloading assemblies...')
db.update_pdbs()
