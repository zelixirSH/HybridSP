
import sys
import argparse
import numpy as np
import multiprocessing as mp


lig_elements = ['C.3', 'C.2', 'C.1', 'C.ar', 'N.3', 'N.2', 'N.1', 'N.ar', 'N.am', 'N.pl3',
                'N.4', 'O.3', 'O.2', 'O.co2', 'S.3', 'S.2', 'P.3', 'C.cat', 'F', 'Cl', 
                'Br', 'I', 'Li', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Cr.th',
                'Cr.oh', 'Mn', 'Fe', 'Co.oh', 'Cu', 'Zn', 'Se', 'Mo', 'Sn']
residues = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU',
            'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
resi_elements = {'ALA': ['N', 'CA', 'O'],      'ARG': ['CD', 'NH1', 'NH2'], 'ASN': ['OD1', 'CG', 'ND2'],
                 'ASP': ['OD1', 'OD2', 'CG'],  'CYS': ['CA', 'CB', 'SG'],   'GLN': ['OE1', 'CG', 'NE2'],
                 'GLU': ['OE1', 'OE2', 'CD'],  'GLY': ['N', 'CA', 'O'],     'HIS': ['CE1', 'CD2', 'NE2'],
                 'ILE': ['CG1', 'CG2', 'CD1'], 'LEU': ['CG', 'CD1', 'CD2'], 'LYS': ['CD', 'CE', 'NZ'],
                 'MET': ['CE', 'CG', 'SD'],    'PHE': ['CD1', 'CD2', 'CZ'], 'PRO': ['CB', 'CD', 'CG'],
                 'SER': ['CB', 'OG', 'N'],     'THR': ['CB', 'OG1', 'CG2'], 'TRP': ['CD2', 'NE1', 'CH2'],
                 'TYR': ['CD1', 'CD2', 'OH'],  'VAL': ['CB', 'CG1', 'CG2']}

def readmol2(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]
    data_sybyl = []
    data_xyz = []
    if not lines:
        print(filename + ' is empty.')
    else:
        atom_index = lines.index('@<TRIPOS>ATOM')
        bond_index = lines.index('@<TRIPOS>BOND')

        for line in lines[atom_index + 1:bond_index]:
            data_sybyl.append(line.split()[5])
            data_xyz.append(line.split()[2:5])

    data_sybyl = np.array(data_sybyl)
    data_xyz = np.array(data_xyz).astype(float)
    return data_sybyl, data_xyz

def readpdb(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        
    data_atom_resi = []
    data_x = []
    data_y = []
    data_z = []
    data_xyz = []
    
    if len(lines) == 0:
        print(filename + ' is empty.')
    else:
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = ''.join(line[12:16].split())
                residue = ''.join(line[17:20].split())
                data_atom_resi.append([atom, residue])
                data_x.append(float(line[30:38]))
                data_y.append(float(line[38:46]))
                data_z.append(float(line[46:54]))
    
    data_atom_resi = np.array(data_atom_resi)
    data_xyz = list(zip(data_x, data_y, data_z))
    data_xyz = np.array(data_xyz)
    
    return data_atom_resi, data_xyz


def get_file_paths(pdb_code, data_path):
    rec_file = str(data_path + '/' + pdb_code + '/' + pdb_code + '_protein.pdb')
    lig_file = str(data_path + '/' + pdb_code + '/' + pdb_code + '_ligand.mol2')
    return rec_file, lig_file

def compute_distances(pos1, pos2):
    return np.linalg.norm(pos1 - pos2)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', type=str, default="data_list", help="Specify the location of data file.")
    parser.add_argument('-dr', type=float, default=0.2, help="The distance interval, default is 0.1.")
    parser.add_argument('-r0', type=float, default=2.0, help="The minimum distance, default is 2.0.")
    parser.add_argument('-rc', type=float, default=6.0, help="The cutoff distance, default is 6.0.")
    parser.add_argument('-weight', action='store_true', help="Enable affinity weighted summation.")
    parser.add_argument('-o', type=str, default="potential_group.npy", help="Output .npy file of the statistical model.")
    parser.add_argument('-path', type=str, default="/data_test/home/lzzheng/wzh/dataset/pdbbind2020/", help="Specify training data path.")

    args = parser.parse_args()
    
    if len(sys.argv) < 2:
        parser.print_help()
        print("\nExample:\npython derive_rdf.py -c lig_select.txt -dr 0.2 -r0 2 -rc 6 -ncpus 50\n")
        sys.exit(0)


    dr = args.dr
    r0 = args.r0
    rc = args.rc
    i0 = int(r0 / dr)
    num_of_mtx = int((rc - r0) / dr)
    data_path = args.path

    with open(args.c, 'r') as rf:
        code_list = []
        pkd_list = []
        for line in rf:
            line = line.strip().split()
            code_list.append(line[0])
            pkd_list.append(float(line[3]))

    pkd_ave = np.mean(pkd_list)
    pkd_tem = [(x - pkd_ave) / pkd_ave for x in pkd_list]

    if args.weight:
        print("Derive Group_W model...")
        p_list = np.exp(pkd_tem)
    else:
        print("Derive Group model...")
        p_list = np.ones_like(pkd_tem)
    
    file_paths = [get_file_paths(pdb_code, data_path) for pdb_code in code_list]
    file_paths_pvalues = [(x[0], x[1], y) for x, y in zip(file_paths, p_list)]

    # atom-3atom
    atom_grp_num_array = np.zeros((23, 20) + (num_of_mtx, num_of_mtx, num_of_mtx)) 

    for rec_file, lig_file, p in file_paths_pvalues:

        rec_atom_resi, rec_xyz = readpdb(rec_file)
        lig_sybyl, lig_xyz = readmol2(lig_file)

        # Compute the distances between sybyl atoms and record them in pair matrices.
        for ind_lig, ele_lig in enumerate(lig_sybyl):
            
            distances_dict = {} # record the distances between ligand and the selected 3 residue atoms
            
            if ele_lig in lig_elements:
                for ind_rec, res_rec in enumerate(rec_atom_resi[:, 1]):
                    if res_rec in resi_elements:
                        ele_rec = rec_atom_resi[ind_rec, 0] # get the element type
                        if ele_rec in resi_elements[res_rec]: # determine whether the element belongs to the seleced ones
                            atom_index = resi_elements[res_rec].index(ele_rec)
                            d = compute_distances(rec_xyz[ind_rec, :], lig_xyz[ind_lig, :])
                            distances_dict[str(atom_index)] = d

                            if len(distances_dict) == 3:
                                if all(r0 < distances_dict[key] < rc for key in ['0', '1', '2']):
                                    outer_row = min(lig_elements.index(ele_lig), 22)
                                    outer_col = residues.index(rec_atom_resi[ind_rec, 1])
                                    inner_row = int((distances_dict['0'] - r0) / dr)
                                    inner_col = int((distances_dict['1'] - r0) / dr)
                                    inner_dep = int((distances_dict['2'] - r0) / dr)

                                    atom_grp_num_array[outer_row, outer_col, inner_row, inner_col, inner_dep] += p

                                distances_dict = {}

    aver_atom_grp_num_array = np.mean(atom_grp_num_array, axis=(0, 1))
    atom_grp_num_array_sum = np.sum(atom_grp_num_array, axis=(0, 1))

    aver_atom_grp_num_array[aver_atom_grp_num_array == 0] = 1
    division_grp = atom_grp_num_array / aver_atom_grp_num_array
    division_grp[division_grp == 0] = 1
    atom_grp_potential_array = - np.log(division_grp)

    np.save(args.o, atom_grp_potential_array)
    print(f"Model has been saved as {args.o}.")

    