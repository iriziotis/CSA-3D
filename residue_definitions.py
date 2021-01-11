from Bio.SeqUtils import IUPACData

AA_3TO1 = IUPACData.protein_letters_3to1_extended.copy()
AA_3TO1[''] = '_'
AA_3TO1['X'] = 'X'

STANDARD_RESIDUES = {
    'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 
    'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
    'MET', 'ASN', 'PRO', 'GLN', 'ARG',
    'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }

RESIDUE_DEFINITIONS = { 
    'ANY.CA', 'ANY.N', 'ANY.O',
    'PTM.CA', 'PTM.C', 'PTM.CB',
    'ALA.C', 'ALA.CA', 'ALA.CB',
    'CYS.CA', 'CYS.CB', 'CYS.SG',
    'ASP.CG', 'ASP.OD1', 'ASP.OD2',
    'GLU.CD', 'GLU.OE1', 'GLU.OE2',
    'PHE.CE1', 'PHE.CE2', 'PHE.CZ',
    'GLY.N', 'GLY.CA', 'GLY.C',
    'HIS.CG', 'HIS.CD2', 'HIS.ND1',
    'ILE.CA', 'ILE.CB', 'ILE.CG1',
    'LYS.CD', 'LYS.CE', 'LYS.NZ',
    'LEU.CA', 'LEU.CB', 'LEU.CG',
    'MET.CG', 'MET.SD', 'MET.CE',
    'ASN.CG', 'ASN.OD1', 'ASN.ND2',
    'PRO.C', 'PRO.CA', 'PRO.O',
    'GLN.CD', 'GLN.OE1', 'GLN.NE2',
    'ARG.CZ', 'ARG.NH1', 'ARG.NH2',
    'SER.CA', 'SER.CB', 'SER.OG',
    'THR.CA', 'THR.CB', 'THR.OG1',
    'TRP.NE1', 'TRP.CZ2', 'TRP.CH2',
    'TYR.CE1', 'TYR.CZ', 'TYR.OH',
    'VAL.CA', 'VAL.CB', 'VAL.CG',
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
        'VAL.CG': 'LEU.CG',
        }
