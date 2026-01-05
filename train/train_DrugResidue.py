
import sys
import argparse
import numpy as np
import multiprocessing as mp


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
                if (residue in resi_elements_dic and atom in resi_elements_dic[residue]) or \
                   (residue in resi_elements_dic['OTH']):
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

def get_contact_num_array(rec_file, lig_file, p=1, dr=0.1, r0=2.0, rc=6.0):

    num_of_mtx = int((rc - r0) / dr)
    local_number_array = [np.zeros((23, 100)) for i in range(num_of_mtx)]
    rec_atom_resi, rec_xyz = readpdb(rec_file)
    lig_sybyl, lig_xyz = readmol2(lig_file)

    # Compute the distances between sybyl atoms and record them in pair matrices.
    for ind_lig, ele_lig in enumerate(lig_sybyl):
        if ele_lig in lig_elements:
            for ind_rec, residue_rec in enumerate(rec_atom_resi[:, 1]):
                if residue_rec in residues:
                    index_rec_part = residues.index(residue_rec)  # get the number of the residue
                    ele_rec = rec_atom_resi[ind_rec, 0] # get the element type
                    if ele_rec in resi_elements[index_rec_part]: # determine whether the element belongs to the seleced ones
                        d = compute_distances(rec_xyz[ind_rec, :], lig_xyz[ind_lig, :])
                        if r0 < d < rc:
                            i = int((d - r0) / dr)
                            row_num = min(lig_elements.index(ele_lig), 22)
                            col_num = 5 * index_rec_part + resi_elements[index_rec_part].index(ele_rec)
                            local_number_array[i][row_num, col_num] += p

                elif residue_rec in resi_elements[20]:
                    d = compute_distances(rec_xyz[ind_rec, :], lig_xyz[ind_lig, :])
                    if r0 < d < rc:
                        i = int((d - r0) / dr)
                        row_num = min(lig_elements.index(ele_lig), 22)
                        col_num = 99
                        local_number_array[i][row_num, col_num] += p
    return local_number_array

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', type=str, default="data_list", help="Specify the location of data file.")
    parser.add_argument('-dr', type=float, default=0.2, help="The distance interval, default is 0.1.")
    parser.add_argument('-r0', type=float, default=2.0, help="The minimum distance, default is 2.0.")
    parser.add_argument('-rc', type=float, default=6.0, help="The cutoff distance, default is 6.0.")
    parser.add_argument('-ncpus', type=int, default=40, help="Number of CPUs, default is 50.")
    parser.add_argument('-weight', action='store_true', help="Enable affinity weighted summation.")
    parser.add_argument('-o', type=str, default="potential_DrugResidue.npy", help="Output .npy file of the statistical model.")
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
    num_processes = args.ncpus
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
        print("Derive DrugResidue_W model...")
        p_list = np.exp(pkd_tem)
    else:
        print("Derive DrugResidue model...")
        p_list = np.ones_like(pkd_tem)
    
    file_paths = [get_file_paths(pdb_code, data_path) for pdb_code in code_list]
    file_paths_pvalues = [(x[0], x[1], y) for x, y in zip(file_paths, p_list)]

    # atom-atom
    number_array = [np.zeros((23, 100)) for i in range(num_of_mtx)] # row: lig; col: rec; depth: distance
    density_array = [np.zeros((23, 100)) for i in range(num_of_mtx)] # row: lig; col: rec; depth: distance
    rdf_array = [np.zeros((23, 100)) for i in range(num_of_mtx)]
    rdf = np.zeros(num_of_mtx)


    with mp.Pool(processes=num_processes) as pool:
        results = pool.starmap(get_contact_num_array, [(rec_file, lig_file, p, dr, r0, rc) for rec_file, lig_file, p in file_paths_pvalues])
    
    number_array = np.sum(results, axis=0)
    number_array_sum = np.sum(number_array, axis=0)

    density_array_sum = number_array_sum / (4/3 * np.pi * rc ** 3 - 4/3 * np.pi * r0 ** 3)

    for i in range(num_of_mtx):
        r = dr * i + r0 + (dr / 2)
        dV = 4 * np.pi * r ** 2 * dr
        density_array[i] = number_array[i] / dV
        rdf_array[i] = np.zeros_like(density_array[i])
        non_zero_mask = density_array_sum != 0
        rdf_array[i][non_zero_mask] = density_array[i][non_zero_mask] / density_array_sum[non_zero_mask]
    
    density_array_sum = np.sum(density_array, axis = 0)
    pair_types_num = int(23 * 100)

    rdf_array_sum = np.sum(rdf_array, axis = 0)
    rdf_array_normalized = [np.zeros((23, 100)) for i in range(num_of_mtx)]
    for i in range(num_of_mtx):
        rdf_array_normalized[i][non_zero_mask] = rdf_array[i][non_zero_mask] / rdf_array_sum[non_zero_mask]
        rdf[i] = np.sum(rdf_array_normalized[i]) / pair_types_num

    ref_rdf = rdf / np.sum(rdf)
    expanded_ref_rdf = ref_rdf[:, np.newaxis, np.newaxis]
    division_atom = rdf_array_normalized / expanded_ref_rdf
    division_atom[division_atom == 0] = 1
    potential_array = - np.log(division_atom)

    
    np.save(args.o, potential_array)
    print(f"Model has been saved as {args.o}.")

