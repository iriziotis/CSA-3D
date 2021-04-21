#!/usr/bin/env python3

import os
import sys
import pickle
from csa3d import Mcsa

def main(mcsa_id, outdir):
    # Initialize database
    db = Mcsa()

    # Get components
    db.get_component_similarities(mcsa_id)


if __name__ == '__main__':
    try:
        outdir = sys.argv[2]
    except IndexError:
        outdir = '.'
    main(int(sys.argv[1]), outdir)
