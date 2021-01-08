import csv
import json
from ast import literal_eval

WORKING_DIR = '/Users/riziotis/ebi/phd'
ASSEMBLIES_DIR = '{}/datasets/assembly_cif'.format(WORKING_DIR)
INFO_JSON = '{}/datasets/catalytic_residues_homologues.json'.format(WORKING_DIR)

uniprot_pdb_mapping = '{}/datasets/sifts/uniprot_pdb.csv'.format(WORKING_DIR)
ligands_json = '{}/datasets/bound_ligands/bound_cognate_similarities/compound_similarities.json'.format(WORKING_DIR)
pdb_ec_csv = '{}/datasets/sifts/pdb_chain_enzyme.csv'.format(WORKING_DIR)

#TODO make me prettier
UNI2PDB = dict()
with open(uniprot_pdb_mapping, 'r') as f:
    next(f)
    for line in csv.reader(f):
        UNI2PDB[line[0]] = line[1].split(';')

#Cognate-bound ligand similarities from PARITY
with open(ligands_json) as infile:
    COMPOUND_SIMILARITIES = {literal_eval(k):float(v) for k, v in json.load(infile).items()}

#PDB ID - E.C. number mapping from SIFTS
with open(pdb_ec_csv) as infile:
    infile.readline()
    reader = csv.reader(infile)
    PDB_EC = {(line[0], line[1]):line[3] for line in reader}
