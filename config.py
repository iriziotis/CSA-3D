import csv
import json
from ast import literal_eval

WORKING_DIR = '/Users/riziotis/ebi/phd'
ASSEMBLIES_DIR = '{}/datasets/assembly_cif'.format(WORKING_DIR)
INFO_JSON = '{}/datasets/catalytic_residues_homologues.json'.format(WORKING_DIR)

uniprot_pdb_mapping_csv = '{}/datasets/sifts/uniprot_pdb.csv'.format(WORKING_DIR)
pdb_ec_mapping_csv = '{}/datasets/sifts/pdb_chain_enzyme.csv'.format(WORKING_DIR)
compound_similarities_json = '{}/datasets/bound_ligands/bound_cognate_similarities/compound_similarities.json'.format(WORKING_DIR)

# UniProt to PDB mapping from sifts
with open(uniprot_pdb_mapping_csv, 'r') as f:
    next(f)
    UNI2PDB = {line[0]: line[1].split(';') for line in csv.reader(f)}

# Cognate-bound ligand similarities from PARITY
with open(compound_similarities_json, 'r') as f:
    COMPOUND_SIMILARITIES = {literal_eval(k): float(v) for k, v in json.load(f).items()}

# PDB ID - E.C. number mapping from SIFTS
with open(pdb_ec_mapping_csv, 'r') as f:
    next(f)
    PDB2EC = {(line[0], line[1]): line[3] for line in csv.reader(f)}
