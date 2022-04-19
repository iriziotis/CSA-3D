from Bio.SeqUtils import IUPACData

AA_1TO3 = IUPACData.protein_letters_1to3_extended.copy()
AA_3TO1 = IUPACData.protein_letters_3to1_extended.copy()
AA_3TO1[''] = '_'
AA_3TO1['X'] = 'X'

STANDARD_RESIDUES = {
    'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 
    'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
    'MET', 'ASN', 'PRO', 'GLN', 'ARG',
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

NUCLEIC = {'DU', 'DT', 'DA', 'DG', 'DC'}

RESIDUE_DEFINITIONS = { 
    'ANY.CA', 'ANY.C', 'ANY.O',
    'PTM.CA', 'PTM.C', 'PTM.CB',
    'GLY.CA', 'GLY.C', 'GLY.O',
    'ALA.CA', 'ALA.C', 'ALA.CB',
    'PRO.CA', 'PRO.C', 'PRO.O',
    'CYS.CA', 'CYS.CB', 'CYS.SG',
    'ASP.CG', 'ASP.OD1', 'ASP.OD2',
    'GLU.CD', 'GLU.OE1', 'GLU.OE2',
    'PHE.CE1', 'PHE.CE2', 'PHE.CZ',
    'HIS.CG', 'HIS.CD2', 'HIS.ND1',
    'ILE.CA', 'ILE.CB', 'ILE.CG1',
    'LYS.CD', 'LYS.CE', 'LYS.NZ',
    'LEU.CA', 'LEU.CB', 'LEU.CG',
    'MET.CG', 'MET.SD', 'MET.CE',
    'ASN.CG', 'ASN.OD1', 'ASN.ND2',
    'GLN.CD', 'GLN.OE1', 'GLN.NE2',
    'ARG.CZ', 'ARG.NH1', 'ARG.NH2',
    'SER.CA', 'SER.CB', 'SER.OG',
    'THR.CA', 'THR.CB', 'THR.OG1',
    'TRP.NE1', 'TRP.CZ2', 'TRP.CH2',
    'TYR.CE1', 'TYR.CZ', 'TYR.OH',
    'VAL.CA', 'VAL.CB', 'VAL.CG1',
    }

EQUIVALENT_RESIDUES = {
        'Glu': 'Asp',
        'Asp': 'Glu',
        'Gln': 'Asn',
        'Asn': 'Gln',
        'Val': {'Leu', 'Ile'},
        'Leu': {'Ile', 'Val'},
        'Ile': {'Val', 'Leu'},
        'Ser': 'Thr',
        'Thr': 'Ser'
        }


EQUIVALENT_ATOMS = {
        'HIS.CD2': 'HIS.ND1',
        'ASP.OD2': 'ASP.OD1',
        'ASN.ND2': 'ASN.OD1',
        'ARG.NH2': 'ARG.NH1',
        'PHE.CE2': 'PHE.CE1',
        'TYR.CE2': 'TYR.CE1',
        'GLU.CD': 'ASP.CG',
        'GLU.OE1': 'ASP.OD1',
        'GLU.OE2': 'ASP.OD1',
        'GLN.CD': 'ASN.CG',
        'GLN.OE1': 'ASN.OD1',
        'GLN.NE2': 'ASN.OD1',
        'SER.CA': 'THR.CA',
        'SER.OG': 'THR.OG1',
        'SER.CB': 'THR.CB',
        'ILE.CA': 'LEU.CA',
        'ILE.CB': 'LEU.CB',
        'ILE.CG1': 'LEU.CG',
        'VAL.CA': 'LEU.CA',
        'VAL.CB': 'LEU.CB',
        'VAL.CG1': 'LEU.CG',
        }

TEMPLATE_FUZZY_RESIDUES = {
        'GLU': 'D',
        'ASP': 'E',
        'GLN': 'N',
        'ASN': 'Q',
        'VAL': 'LI',
        'LEU': 'VI',
        'ILE': 'LV',
        'SER': 'TY',
        'THR': 'SY',
        'TYR': 'ST '
        }

TEMPLATE_FUZZY_ATOMS = { 
    # TODO Need to cross check these, also see non-fuzzy case
    'ANY': {'CA': [100, ''], 'C': [100, ''], 'O': [100, '']},
    'PTM': {'CA': [100, ''], 'C': [100, ''],  'CB': [100, '']},
    'ALA': {'C' : [0, ''], 'CA' : [0, ''], 'CB' : [0, '']},
    'CYS': {'CA' : [0, ''], 'CB' : [0, ''], 'SG' : [0, '']},
    'ASP': {'CG': [3, 'E'], 'OD1': [3, 'E'], 'OD2': [3, 'E']},
    'GLU': {'CD': [3, 'D'], 'OE1': [3, 'D'], 'OE2': [3, 'D']},
    'PHE': {'CE1': [3, ''], 'CE2': [3, ''], 'CZ': [0, '']},
    'GLY': {'O' : [0, ''], 'CA' : [0, ''], 'C' : [0, '']},
    'HIS': {'CG': [0, ''], 'CD2':[8, ''], 'ND1': [8, '']},
    'ILE': {'CA': [0, 'VL'], 'CB': [0, 'VL'], 'CG1': [3, 'VL']},
    'LYS': {'CD' : [0, ''], 'CE' : [0, ''], 'NZ' : [0, '']},
    'LEU': {'CA': [0, 'VI'], 'CB': [0, 'VI'], 'CG': [3, 'VI']},
    'MET': {'CG' : [0, ''], 'SD' : [0, ''], 'CE' : [0, '']},
    'ASN': {'CG': [3, 'Q'], 'OD1': [1, 'Q'] , 'ND2': [1, 'Q']},
    'GLN': {'CD': [3, 'N'], 'OE1': [1, 'N'], 'NE2': [1, 'N']},
    'PRO': {'C' : [0, ''], 'CA' : [0, ''], 'O' : [0, '']},
    'ARG': {'CZ': [0, ''], 'NH1': [3, ''], 'NH2': [3, '']},
    'SER': {'CA': [3, 'TY'], 'CB':[3, 'TY'], 'OG': [1, 'TY']},
    'THR': {'CA': [3, 'SY'], 'CB': [3, 'SY'], 'OG1': [1, 'SY']},
    'TYR': {'CE1': [3, 'TS'], 'CZ': [3, 'TS'], 'OH': [1, 'TS']},
    'TRP': {'NE1' : [0, ''], 'CZ2' : [0, ''], 'CH2' : [0, '']},
    'VAL': {'CA': [0, 'LI'], 'CB': [0, 'LI'], 'CG': [3, 'LI']},
    }
