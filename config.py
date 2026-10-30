"""
Configuration for CSA-3D data paths.

Uses flexible paths from csa3d.config, falling back to hardcoded paths for legacy compatibility.
"""

import os
import csv
import json
from ast import literal_eval
from collections import defaultdict

# Try to use csa3d's flexible config first
try:
    from csa3d.config import get_data_dir
    _data_dir = get_data_dir()
except (ImportError, Exception):
    # Fallback to hardcoded paths for legacy/offline use
    _data_dir = '/Users/riziotis/ebi/phd/data/'
    if not os.path.isdir(_data_dir):
        _data_dir = '/nfs/research/thornton/riziotis/research/phd/data/'

# Build all paths relative to the data directory
WORKING_DIR = _data_dir
ASSEMBLIES_DIR = os.path.join(WORKING_DIR, 'pdbe/assembly_cif/cif')
KEGG_MOLS_DIR = os.path.join(WORKING_DIR, 'kegg/sdf')
CHEBI_MOLS_DIR = os.path.join(WORKING_DIR, 'chebi/sdf')
HET_MOLS_DIR = os.path.join(WORKING_DIR, 'pdbe/het_sdf')
CAT_RES_INFO = os.path.join(WORKING_DIR, 'mcsa/catalytic_residues_homologues.json')
REF_RES_INFO = os.path.join(WORKING_DIR, 'mcsa/catalytic_residues_references.json')
MCSA_ENTRY_INFO = os.path.join(WORKING_DIR, 'mcsa/entry_info/mcsa_entry_info.*.json')

# SIFTS
uniprot_pdb_mapping_csv = os.path.join(WORKING_DIR, 'sifts/uniprot_pdb.csv')
pdb_uniprot_mapping_csv = os.path.join(WORKING_DIR, 'sifts/pdb_chain_uniprot.csv')
pdb_ec_mapping_csv = os.path.join(WORKING_DIR, 'sifts/pdb_chain_enzyme.csv')
# KEGG
kegg_reaction_mapping_csv = os.path.join(WORKING_DIR, 'kegg/kegg_reaction_components.csv')
# RHEA
rhea_reaction_mapping_csv = os.path.join(WORKING_DIR, 'rhea/rhea_reaction_components.csv')
# Ligands
pdbe_ions_csv = os.path.join(WORKING_DIR, 'pdbe/pdb_ions.csv')
pdbe_cofactors_csv = os.path.join(WORKING_DIR, 'pdbe/pdb_cofactors.csv')
mcsa_cofactors_csv = os.path.join(WORKING_DIR, 'mcsa/mcsa_cofactors.csv')
crystallization_hets_csv = os.path.join(WORKING_DIR, 'pdbe/crystallization_hets.csv')
# Resolution
pdb_resolution = os.path.join(WORKING_DIR, 'pdbe/resolu.idx')

# Load data files with error handling
def _safe_load_csv(filepath, dict_key=None, skip_header=False):
    """Safely load CSV file with informative error message."""
    try:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Data file not found: {filepath}")

        with open(filepath, 'r') as f:
            if skip_header:
                next(f)
            if dict_key is not None and skip_header:
                next(f)  # Skip second header if needed
            reader = csv.reader(f)
            return reader if dict_key is None else {line[0]: line[dict_key] for line in reader if line}
        return {}
    except Exception as e:
        raise FileNotFoundError(
            f"Could not load {filepath}. "
            f"Make sure data is set up with: csa3d.setup_data()\n"
            f"Error: {e}"
        )

# UniProt to PDB mapping from sifts
try:
    with open(uniprot_pdb_mapping_csv, 'r') as f:
        next(f)
        UNI2PDB = {line[0]: line[1].split(';') for line in csv.reader(f) if len(line) >= 2}
except Exception as e:
    print(f"Warning: Could not load {uniprot_pdb_mapping_csv}: {e}")
    UNI2PDB = {}

# PDB to UniProt mapping from sifts
try:
    with open(pdb_uniprot_mapping_csv, 'r') as f:
        next(f)
        next(f)
        PDB2UNI = {(line[0], line[1]): line[2] for line in csv.reader(f) if len(line) >= 3}
except Exception as e:
    print(f"Warning: Could not load {pdb_uniprot_mapping_csv}: {e}")
    PDB2UNI = {}

# PDB ID - E.C. number mapping from SIFTS
try:
    with open(pdb_ec_mapping_csv, 'r') as f:
        next(f)
        PDB2EC = {(line[0], line[1]): line[3] for line in csv.reader(f) if len(line) >= 4}
except Exception as e:
    print(f"Warning: Could not load {pdb_ec_mapping_csv}: {e}")
    PDB2EC = {}

# KEGG ec - reaction components mapping
try:
    with open(kegg_reaction_mapping_csv, 'r') as f:
        KEGG_EC_REACTION = {line[0]: [line[1].split(';')] for line in csv.reader(f) if len(line) >= 2}
except Exception as e:
    print(f"Warning: Could not load {kegg_reaction_mapping_csv}: {e}")
    KEGG_EC_REACTION = {}

# RHEA ec - reaction components mapping
try:
    with open(rhea_reaction_mapping_csv, 'r') as f:
        next(f)
        RHEA_EC_REACTION = {line[2]: line[3].split(';') for line in csv.reader(f) if len(line) >= 4}
except Exception as e:
    print(f"Warning: Could not load {rhea_reaction_mapping_csv}: {e}")
    RHEA_EC_REACTION = {}

# Crystallization artefact HETs set
try:
    with open(crystallization_hets_csv, 'r') as f:
        CRYSTALLIZATION_HETS = {line.strip() for line in f if line.strip()}
except Exception as e:
    print(f"Warning: Could not load {crystallization_hets_csv}: {e}")
    CRYSTALLIZATION_HETS = set()

# Ions sets
METALS = set()
REACTIVE_NONMETALS = set()
NOBLE_GASES = set()
try:
    with open(pdbe_ions_csv, 'r') as f:
        next(f)
        for line in csv.reader(f):
            if len(line) >= 4:
                if line[3] == 'METAL':
                    METALS.add(line[0])
                if line[3] == 'REACTIVE_NONMETAL':
                    REACTIVE_NONMETALS.add(line[0])
                if line[3] == 'NOBLE_GAS':
                    NOBLE_GASES.add(line[0])
except Exception as e:
    print(f"Warning: Could not load {pdbe_ions_csv}: {e}")

# PDB ID - co-factors mapping
PDBID_COFACTORS = defaultdict(set)
try:
    with open(pdbe_cofactors_csv, 'r') as f:
        next(f)
        for line in csv.reader(f, quotechar='"'):
            if len(line) >= 5:
                PDBID_COFACTORS[line[0]].add(line[4])
except Exception as e:
    print(f"Warning: Could not load {pdbe_cofactors_csv}: {e}")

# MCSA EC - co-factors mapping
MCSA_EC_COFACTORS = defaultdict(set)
try:
    with open(mcsa_cofactors_csv, 'r') as f:
        next(f)
        for line in csv.reader(f, quotechar='"'):
            if len(line) >= 5:
                for ec in line[4].split(';'):
                    MCSA_EC_COFACTORS[ec].add(line[1])
except Exception as e:
    print(f"Warning: Could not load {mcsa_cofactors_csv}: {e}")

# PDB_ID - resolution
try:
    with open(pdb_resolution, 'r') as f:
        PDBID_RESOLUTION = {line[0].strip().lower(): line[1].strip() for line in csv.reader(f, delimiter=';') if len(line) >= 2}
except Exception as e:
    print(f"Warning: Could not load {pdb_resolution}: {e}")
    PDBID_RESOLUTION = {}


# Crystallization artefact HETs set
try:
    with open(crystallization_hets_csv, 'r') as f:
        CRYSTALLIZATION_HETS = {line.strip() for line in f}
except FileNotFoundError as e:
    print(f"Warning: {e}")
    CRYSTALLIZATION_HETS = set()

# Ions sets
METALS = set()
REACTIVE_NONMETALS = set()
NOBLE_GASES = set()
try:
    with open(pdbe_ions_csv, 'r') as f:
        next(f)
        for line in csv.reader(f):
            if line[3] == 'METAL':
                METALS.add(line[0])
            if line[3] == 'REACTIVE_NONMETAL':
                REACTIVE_NONMETALS.add(line[0])
            if line[3] == 'NOBLE_GAS':
                NOBLE_GASES.add(line[0])
except FileNotFoundError as e:
    print(f"Warning: {e}")

# PDB ID - co-factors mapping
PDBID_COFACTORS = defaultdict(set)
try:
    with open(pdbe_cofactors_csv, 'r') as f:
        next(f)
        for line in csv.reader(f, quotechar='"'):
            PDBID_COFACTORS[line[0]].add(line[4])
except FileNotFoundError as e:
    print(f"Warning: {e}")

# MCSA EC - co-factors mapping
MCSA_EC_COFACTORS = defaultdict(set)
try:
    with open(mcsa_cofactors_csv, 'r') as f:
        next(f)
        for line in csv.reader(f, quotechar='"'):
            for ec in line[4].split(';'):
                MCSA_EC_COFACTORS[ec].add(line[1])
except FileNotFoundError as e:
    print(f"Warning: {e}")

# PDB_ID - resolution
try:
    with open(pdb_resolution, 'r') as f:
        PDBID_RESOLUTION = {line[0].strip().lower(): line[1].strip() for line in csv.reader(f, delimiter=';')}
except FileNotFoundError as e:
    print(f"Warning: {e}")
    PDBID_RESOLUTION = {}
