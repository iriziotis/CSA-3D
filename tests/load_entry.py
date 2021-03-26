#!/usr/bin/env python3

import sys
import os
import pickle
import Bio.PDB

def main(mcsa_id):
    with open('entries/csa3d_{}.ent'.format(str(mcsa_id).zfill(4)), 'rb') as f:
        entry = pickle.load(f)

    outdir = './out/mcsa_{}'.format(str(mcsa_id).zfill(4))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for pdbsite in entry.pdbsites:
        rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, transform=True)
        per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite, isolate_residue=True)
        print(pdbsite.id, per_res_rms)
        pdbsite.write_pdb(outdir=outdir, write_hets=True, func_atoms_only=False)

    # To test for wrong superpositions (entry 925)
    #a = entry.get_pdbsite('4pw9_A-A-A')
    #b = entry.get_pdbsite('3hbq_A-A-A')
    #
    #i = None
    #rot, tran, rms, rms_all = a.fit(b, exclude=i)
    #print('Excluding', i, rms, rms_all)
    #b.structure.transform(rot, tran)
    #
    #a.write_pdb(outfile='a.pdb', write_hets=False, func_atoms_only=True)
    #b.write_pdb(outfile='b.pdb', write_hets=False, func_atoms_only=True)
    #
    #print(a.per_residue_rms(b, isolate_residue=True))

if __name__ == '__main__':
    main(int(sys.argv[1]))
                
