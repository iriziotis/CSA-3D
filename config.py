import os
import csv
import json
from ast import literal_eval
from collections import defaultdict

WORKING_DIR = '/Users/riziotis/ebi/phd'
if not os.path.isdir(WORKING_DIR):
    WORKING_DIR = '/nfs/research1/thornton/riziotis/research/phd'
ASSEMBLIES_DIR = '{}/datasets/assembly_cif'.format(WORKING_DIR)
CAT_RES_INFO = '{}/datasets/catalytic_residues_homologues.json'.format(WORKING_DIR)
MCSA_ENTRY_INFO = '{}/datasets/mcsa_entry_info/mcsa_entry_info.*.json'.format(WORKING_DIR)

uniprot_pdb_mapping_csv = '{}/datasets/sifts/uniprot_pdb.csv'.format(WORKING_DIR)
pdb_uniprot_mapping_csv = '{}/datasets/sifts/pdb_chain_uniprot.csv'.format(WORKING_DIR)
pdb_ec_mapping_csv = '{}/datasets/sifts/pdb_chain_enzyme.csv'.format(WORKING_DIR)
compound_similarities_json = '{}/datasets/bound_ligands/parity_data/compound_similarities.json'.format(WORKING_DIR)
pdbe_cofactors_csv = '{}/datasets/bound_ligands/pdbe-kb_data/pdb_cofactors.csv'.format(WORKING_DIR)
mcsa_cofactors_csv = '{}/datasets/bound_ligands/mcsa_data/mcsa_cofactors.csv'.format(WORKING_DIR)
crystallization_hets = '{}/datasets/crystallization_hets.csv'.format(WORKING_DIR)

# UniProt to PDB mapping from sifts
with open(uniprot_pdb_mapping_csv, 'r') as f:
    next(f)
    UNI2PDB = {line[0]: line[1].split(';') for line in csv.reader(f)}

# PDB to UniProt mapping from sifts
with open(pdb_uniprot_mapping_csv, 'r') as f:
    next(f)
    next(f)
    PDB2UNI = {(line[0], line[1]): line[2] for line in csv.reader(f)}

# Cognate-bound ligand similarities from PARITY
with open(compound_similarities_json, 'r') as f:
    COMPOUND_SIMILARITIES = {literal_eval(k): float(v) for k, v in json.load(f).items()}

# PDB ID - E.C. number mapping from SIFTS
with open(pdb_ec_mapping_csv, 'r') as f:
    next(f)
    PDB2EC = {(line[0], line[1]): line[3] for line in csv.reader(f)}

# PDB ID - co-factors mapping
PDBID_COFACTORS = defaultdict(set)
with open(pdbe_cofactors_csv, 'r') as f:
    next(f)
    for line in csv.reader(f, quotechar='"'):
        PDBID_COFACTORS[line[0]].add(line[4])

# Metal co-factor set
with open(mcsa_cofactors_csv, 'r') as f:
    next(f)
    METAL_COFACTORS = {line[1] for line in csv.reader(f, quotechar='"') if line[3]=="True"}

# Crystallization artefact HETs set
with open(crystallization_hets, 'r') as f:
    CRYSTALLIZATION_HETS = {line.strip() for line in f}

