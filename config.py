import os
import csv
import json
from ast import literal_eval
from collections import defaultdict

WORKING_DIR = '/Users/riziotis/ebi/phd/datasets/'
if not os.path.isdir(WORKING_DIR):
    WORKING_DIR = '/nfs/research1/thornton/riziotis/research/phd/datasets/'
ASSEMBLIES_DIR = WORKING_DIR + 'pdbe/assembly_cif'
KEGG_MOLS_DIR = WORKING_DIR + 'kegg/mol'
CHEBI_MOLS_DIR = WORKING_DIR + 'chebi/sdf'
HET_MOLS_DIR = WORKING_DIR + 'pdbe/het_sdf'
CAT_RES_INFO = WORKING_DIR + 'mcsa/catalytic_residues_homologues.json'
MCSA_ENTRY_INFO = WORKING_DIR + 'mcsa/entry_info/mcsa_entry_info.*.json'

# SIFTS
uniprot_pdb_mapping_csv = WORKING_DIR + 'sifts/uniprot_pdb.csv'
pdb_uniprot_mapping_csv = WORKING_DIR + 'sifts/pdb_chain_uniprot.csv'
pdb_ec_mapping_csv = WORKING_DIR + 'sifts/pdb_chain_enzyme.csv'
# KEGG
kegg_reaction_mapping_csv = WORKING_DIR + 'kegg/kegg_reaction_components.csv'
# RHEA
rhea_reaction_mapping_csv = WORKING_DIR + 'rhea/rhea_reaction_components.csv'
# Ligands
pdbe_ions_csv = WORKING_DIR + 'pdbe/pdb_ions.csv'
pdbe_cofactors_csv = WORKING_DIR + 'pdbe/pdb_cofactors.csv'
mcsa_cofactors_csv = WORKING_DIR + 'mcsa/mcsa_cofactors.csv'
crystallization_hets_csv = WORKING_DIR + 'pdbe/crystallization_hets.csv'

# UniProt to PDB mapping from sifts
with open(uniprot_pdb_mapping_csv, 'r') as f:
    next(f)
    UNI2PDB = {line[0]: line[1].split(';') for line in csv.reader(f)}

# PDB to UniProt mapping from sifts
with open(pdb_uniprot_mapping_csv, 'r') as f:
    next(f)
    next(f)
    PDB2UNI = {(line[0], line[1]): line[2] for line in csv.reader(f)}

# PDB ID - E.C. number mapping from SIFTS
with open(pdb_ec_mapping_csv, 'r') as f:
    next(f)
    PDB2EC = {(line[0], line[1]): line[3] for line in csv.reader(f)}

# KEGG ec - reaction components mapping
with open(kegg_reaction_mapping_csv, 'r') as f:
    KEGG_EC_REACTION = {line[0]: [line[1].split(';')] for line in csv.reader(f)}

# RHEA ec - reaction components mapping
with open(rhea_reaction_mapping_csv, 'r') as f:
    next(f)
    RHEA_EC_REACTION = {line[2]: line[3].split(';') for line in csv.reader(f)}

# Crystallization artefact HETs set
with open(crystallization_hets_csv, 'r') as f:
    CRYSTALLIZATION_HETS = {line.strip() for line in f}

# Ions sets
METALS = set()
REACTIVE_NONMETALS = set()
NOBLE_GASES = set()
with open(pdbe_ions_csv, 'r') as f:
    next(f)
    for line in csv.reader(f):
        if line[3] == 'METAL':
            METALS.add(line[0])
        if line[3] == 'REACTIVE_NONMETAL':
            REACTIVE_NONMETALS.add(line[0])
        if line[3] == 'NOBLE_GAS':
            NOBLE_GASES.add(line[0])

# PDB ID - co-factors mapping
PDBID_COFACTORS = defaultdict(set)
with open(pdbe_cofactors_csv, 'r') as f:
    next(f)
    for line in csv.reader(f, quotechar='"'):
        PDBID_COFACTORS[line[0]].add(line[4])

# MCSA EC - co-factors mapping
MCSA_EC_COFACTORS = defaultdict(set)
with open(mcsa_cofactors_csv, 'r') as f:
    next(f)
    for line in csv.reader(f, quotechar='"'):
        for ec in line[4].split(';'):
            MCSA_EC_COFACTORS[ec].add(line[1])

