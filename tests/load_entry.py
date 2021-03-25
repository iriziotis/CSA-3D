#!/usr/bin/env python3

import sys
import os
import pickle
import Bio.PDB

def main(mcsa_id):
    with open('csa3d_{}.ent'.format(str(mcsa_id).zfill(4)), 'rb') as f:
        entry = pickle.load(f)

    outdir = './out/mcsa_{}'.format(str(mcsa_id).zfill(4))
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #a = entry.get_pdbsite('1q52_B-B-B-B-B-B')
    #b = entry.get_pdbsite('4q1i_C-C-C-C-C-C')
    #rot, tran, rms, rms_all = a.fit(b)
    #print(rms, rms_all)
    #for i, _ in enumerate(a.residues):
    #    rot, tran, rms, rms_all = a.fit(b, exclude=i)
    #    b.structure.transform(rot, tran)
    #    b.write_pdb(outfile=f'{i}.pdb', write_hets=False, func_atoms_only=True)
    #b.write_pdb(outfile='1q52site.pdb', write_hets=False, func_atoms_only=False)
    
    for pdbsite in entry.pdbsites:
        rot, tran, rms, rms_all = pdbsite.reference_site.fit(pdbsite, transform=True)
        per_res_rms = pdbsite.reference_site.per_residue_rms(pdbsite, isolate_residue=True)
        print(pdbsite.id, per_res_rms)
        pdbsite.write_pdb(outdir=outdir, write_hets=True, func_atoms_only=False)


    
if __name__ == '__main__':
    main(int(sys.argv[1]))
                
