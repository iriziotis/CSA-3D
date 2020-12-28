import csv

WORKING_DIR = '/Users/riziotis/ebi/phd'
INFO_JSON = '{}/datasets/catalytic_residues_homologues.json'.format(WORKING_DIR)
uniprot_pdb_mapping = '{}/datasets/sifts/uniprot_pdb.csv'.format(WORKING_DIR)
UNI2PDB = dict()
with open(uniprot_pdb_mapping, 'r') as f:
    next(f)
    for line in csv.reader(f):
        UNI2PDB[line[0]] = line[1].split(';')
