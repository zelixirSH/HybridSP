from scripts.load_file import readmol2, readpdb
from scripts.util import lig_elements, residues_grp, resi_elements_grp
import numpy as np


def compute_distances(pos1, pos2):
    return np.linalg.norm(pos1 - pos2)


def get_atom_grp_num_array(rec_file, lig_file, dr=0.2, r0=2.0, rc=6.0):
    num_of_mtx = int((rc - r0) / dr)
    outer_shape = (23, 20)
    inner_shape = (num_of_mtx, num_of_mtx, num_of_mtx)

    lig_mols_list, min_xyz, max_xyz = readmol2(lig_file)
    rec_atom_resi, rec_xyz = readpdb(rec_file, min_xyz, max_xyz)
    
    local_grp_number_arrays_list = []
    for lig_name, lig_sybyl, lig_xyz in lig_mols_list:
        local_grp_number_array = np.zeros(outer_shape + inner_shape)
        
        for ind_lig, ele_lig in enumerate(lig_sybyl):
            distances_dict = {}  # record the distances between ligand and the selected 3 residue atoms
            
            if ele_lig in lig_elements:
                for ind_rec, residue_type in enumerate(rec_atom_resi[:, 1]):
                    if residue_type in residues_grp:
                        ele_rec = rec_atom_resi[ind_rec, 0]
                        
                        if ele_rec in resi_elements_grp[residue_type]:
                            atom_index = resi_elements_grp[residue_type].index(ele_rec)
                            d = compute_distances(rec_xyz[ind_rec, :], lig_xyz[ind_lig, :])
                            distances_dict[str(atom_index)] = d

                            if len(distances_dict) == 3:
                                if all(r0 < distances_dict[key] < rc for key in ['0', '1', '2']):
                                    outer_row = min(lig_elements.index(ele_lig), 22)
                                    outer_col = residues_grp.index(residue_type)
                                    inner_row = int((distances_dict['0'] - r0) / dr)
                                    inner_col = int((distances_dict['1'] - r0) / dr)
                                    inner_dep = int((distances_dict['2'] - r0) / dr)
                                    
                                    local_grp_number_array[outer_row, outer_col, inner_row, inner_col, inner_dep] += 1
                                distances_dict = {}
        local_grp_number_arrays_list.append(local_grp_number_array)
    
    return local_grp_number_arrays_list


def group_score(rec_file, lig_file, atom_grp_potential_array, dr=0.2, r0=2.0, rc=6.0):
    grp_number_array_list = get_atom_grp_num_array(rec_file, lig_file, dr, r0, rc)
    
    score_list = []
    for n in range(len(grp_number_array_list)):
        grp_number_array = grp_number_array_list[n]
        atom_energy_array = grp_number_array * atom_grp_potential_array
        potential = np.sum(atom_energy_array)
        score_list.append(potential)
    return np.array(score_list, dtype=float)