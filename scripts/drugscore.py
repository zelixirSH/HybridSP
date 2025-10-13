from scripts.util import drugscore_elements
from scripts.load_file import read_lig_mol2_drugscore, read_rec_mol2_drugscore
import subprocess
import numpy as np
import uuid
import os


def compute_distances(pos1, pos2):
    return np.linalg.norm(pos1 - pos2)


def get_contact_num_array(rec_file, lig_file, dr=0.1, r0=2.0, rc=6.0):
    num_of_mtx = int((rc - r0) / dr)
    lig_mols_list, min_xyz, max_xyz = read_lig_mol2_drugscore(lig_file)
    rec_sybyl, rec_xyz = read_rec_mol2_drugscore(rec_file)
    local_number_arrays_list = []
    lig_name_list = []
    for lig_name, lig_sybyl, lig_xyz in lig_mols_list:
        local_number_array = [np.zeros((25, 25)) for i in range(num_of_mtx)]
        lig_name_list.append(lig_name)
        # Compute the distances between sybyl atoms and record them in pair matrices.
        for ind_lig, ele_lig in enumerate(lig_sybyl):
            for ind_rec, ele_rec in enumerate(rec_sybyl):
                if ele_rec in drugscore_elements:
                    row_num = drugscore_elements.index(ele_lig)
                    col_num = drugscore_elements.index(ele_rec)
                    if row_num >= 24:
                        row_num = 24
                    if col_num >= 24:
                        col_num = 24

                    d = compute_distances(rec_xyz[ind_rec, :], lig_xyz[ind_lig, :])
                    
                    if r0 < d < rc:
                        i = int((d-r0) / dr)
                        local_number_array[i][row_num, col_num] += 1

        local_number_arrays_list.append(local_number_array)
    return lig_name_list, local_number_arrays_list



def drugscore_score(rec_file, lig_file, potential_array, tmp='./tmp', dr=0.1, r0=2.0, rc=6.0):
    os.makedirs(tmp, exist_ok=True)
    tmp_pro_mol2 = f'{tmp}/{uuid.uuid4().hex[:20]}.mol2'
    cmd_obabel = f'obabel {rec_file} -O {tmp_pro_mol2}'

    try:
        obabel_result = subprocess.run(cmd_obabel, shell=True, capture_output=True, 
                                       text=True, check=True)
    
    except subprocess.CalledProcessError as e:
        print(f"DrugScore error while converting {pro_pdb} to mol2: {e}")
        return None
    
    lig_name_list, atom_number_arrays_list = get_contact_num_array(tmp_pro_mol2, lig_file, dr, r0, rc)
    
    score_list = []
    for n, lig_name in enumerate(lig_name_list):
        atom_number_array = atom_number_arrays_list[n]
        atom_energy_array = atom_number_array * potential_array
        potential = np.sum(atom_energy_array)
        score_list.append(potential)
    
    os.remove(tmp_pro_mol2)

    return lig_name_list, np.array(score_list, dtype=float)