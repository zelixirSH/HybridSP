lig_elements = ['C.3', 'C.2', 'C.1', 'C.ar', 'N.3', 'N.2', 'N.1', 'N.ar', 'N.am', 'N.pl3',
                'N.4', 'O.3', 'O.2', 'O.co2', 'S.3', 'S.2', 'P.3', 'C.cat', 'F', 'Cl',
                'Br', 'I', 'Li', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Cr.th',
                'Cr.oh', 'Mn', 'Fe', 'Co.oh', 'Cu', 'Zn', 'Se', 'Mo', 'Sn']


residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 
            'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'GLY']


resi_elements = [['N', 'CA', 'CB', 'O', 'C'],      ['N', 'O', 'CD', 'NH1', 'NH2'], 
                 ['N', 'O', 'OD1', 'CG', 'ND2'],   ['N', 'O', 'OD1', 'OD2', 'CG'],
                 ['N', 'CA', 'O', 'CB', 'SG'],     ['N', 'O', 'OE1', 'CG', 'NE2'],
                 ['N', 'O', 'OE1', 'OE2', 'CG'],   ['N', 'O', 'CB', 'ND1', 'NE2'],
                 ['N', 'O', 'CG1', 'CG2', 'CD1'],  ['N', 'O', 'CG', 'CD1', 'CD2'],
                 ['N', 'O', 'CD', 'CE', 'NZ'],     ['N', 'O', 'CE', 'CG', 'SD'],
                 ['N', 'O', 'CD1', 'CD2', 'CZ'],   ['N', 'O', 'CB', 'CD', 'CG'],
                 ['N', 'CA', 'CB', 'OG', 'O'],     ['N', 'O', 'CB', 'OG1', 'CG2'],
                 ['N', 'O', 'CD2', 'NE1', 'CH2'],  ['N', 'O', 'CD1', 'CD2', 'OH'],
                 ['N', 'O', 'CB', 'CG1', 'CG2'],   ['N', 'CA', 'O', 'C'],
                 ['F', 'CL', 'BR', 'I', 'LI', 'NA', 'MG', 'AL', 'SI', 'K', 'CA', 'CR', 
                  'MN', 'FE', 'CO', 'CU', 'ZN', 'SE', 'MO', 'SN']]


resi_elements_dic = {'ALA': ['N', 'CA', 'CB', 'O', 'C'],      'ARG': ['N', 'O', 'CD', 'NH1', 'NH2'], 
                     'ASN': ['N', 'O', 'OD1', 'CG', 'ND2'],   'ASP': ['N', 'O', 'OD1', 'OD2', 'CG'], 
                     'CYS': ['N', 'CA', 'O', 'CB', 'SG'],     'GLN': ['N', 'O', 'OE1', 'CG', 'NE2'], 
                     'GLU': ['N', 'O', 'OE1', 'OE2', 'CG', 'CD'],   'HIS': ['N', 'O', 'CB', 'ND1', 'NE2', 'CE1', 'CD2'],
                     'ILE': ['N', 'O', 'CG1', 'CG2', 'CD1'],  'LEU': ['N', 'O', 'CG', 'CD1', 'CD2'], 
                     'LYS': ['N', 'O', 'CD', 'CE', 'NZ'],     'MET': ['N', 'O', 'CE', 'CG', 'SD'],
                     'PHE': ['N', 'O', 'CD1', 'CD2', 'CZ'],   'PRO': ['N', 'O', 'CB', 'CD', 'CG'], 
                     'SER': ['N', 'CA', 'CB', 'OG', 'O'],     'THR': ['N', 'O', 'CB', 'OG1', 'CG2'], 
                     'TRP': ['N', 'O', 'CD2', 'NE1', 'CH2'],  'TYR': ['N', 'O', 'CD1', 'CD2', 'OH'], 
                     'VAL': ['N', 'O', 'CB', 'CG1', 'CG2'],   'GLY': ['N', 'CA', 'O', 'C'],
                     'OTH': ['F', 'CL', 'BR', 'I', 'LI', 'NA', 'MG', 'AL', 'SI', 'K', 'CA', 'CR', 
                             'MN', 'FE', 'CO', 'CU', 'ZN', 'SE', 'MO', 'SN']}


residues_grp = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


resi_elements_grp = {
    'ALA': ['N', 'CA', 'O'], 'ARG': ['CD', 'NH1', 'NH2'], 'ASN': ['OD1', 'CG', 'ND2'],
    'ASP': ['OD1', 'OD2', 'CG'], 'CYS': ['CA', 'CB', 'SG'], 'GLN': ['OE1', 'CG', 'NE2'],
    'GLU': ['OE1', 'OE2', 'CD'], 'GLY': ['N', 'CA', 'O'], 'HIS': ['CE1', 'CD2', 'NE2'],
    'ILE': ['CG1', 'CG2', 'CD1'], 'LEU': ['CG', 'CD1', 'CD2'], 'LYS': ['CD', 'CE', 'NZ'],
    'MET': ['CE', 'CG', 'SD'], 'PHE': ['CD1', 'CD2', 'CZ'], 'PRO': ['CB', 'CD', 'CG'],
    'SER': ['CB', 'OG', 'N'], 'THR': ['CB', 'OG1', 'CG2'], 'TRP': ['CD2', 'NE1', 'CH2'],
    'TYR': ['CD1', 'CD2', 'OH'], 'VAL': ['CB', 'CG1', 'CG2']
    }


drugscore_elements = ['C.3', 'C.2', 'C.1', 'C.ar', 'N.3', 'N.2', 'N.1', 'N.ar', 'N.am', 'N.pl3',
                      'N.4', 'O.3', 'O.2', 'O.co2', 'S.3', 'S.2', 'S.O', 'P.3', 'C.cat', 'S.O2', 
                      'F', 'Cl', 'Br', 'I', 'Li', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Cr.th',
                      'Cr.oh', 'Mn', 'Fe', 'Co.oh', 'Cu', 'Zn', 'Se', 'Mo', 'Sn']