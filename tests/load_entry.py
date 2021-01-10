#!/usr/bin/env python3
import sys
sys.path.append('/Users/riziotis/ebi/phd/src/')
import pickle

def main(mcsa_id):
    with open('entry_{}.ent'.format(mcsa_id), 'rb') as f:
        entry = pickle.load(f)

    print('Fitting')
    for pdbsite in entry.pdbsites:
        print(pdbsite.id, pdbsite)
        if pdbsite.is_reference:
            pdbsite.write_pdb(outdir='./out', write_hets=True)
            continue
        try:
            pdbsite.reference_site.fit(pdbsite, transform=True)
        except Exception as e:
            print(e)
            continue
        pdbsite.write_pdb(outdir='./out', write_hets=True)


if __name__ == '__main__':
    main(int(sys.argv[1]))

