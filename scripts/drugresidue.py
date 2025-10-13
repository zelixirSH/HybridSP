from scripts.util import residues, lig_elements, resi_elements
from scripts.load_file import readmol2, readpdb
import numpy as np


def compute_distances(pos1, pos2):
    return np.linalg.norm(pos1 - pos2)


def get_contact_num_array(rec_file, lig_file, dr=0.2, r0=2.0, rc=6.0):

    num_of_mtx = int((rc - r0) / dr)
    lig_mols_list, min_xyz, max_xyz = readmol2(lig_file)
    rec_atom_resi, rec_xyz = readpdb(rec_file, min_xyz, max_xyz)
    local_number_arrays_list = []
    lig_name_list = []
    for lig_name, lig_sybyl, lig_xyz in lig_mols_list:
        local_number_array = [np.zeros((23, 100)) for i in range(num_of_mtx)]
        lig_name_list.append(lig_name)
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
                                local_number_array[i][row_num, col_num] += 1

                    elif residue_rec in resi_elements[20]:
                        d = compute_distances(rec_xyz[ind_rec, :], lig_xyz[ind_lig, :])
                        if r0 < d < rc:
                            i = int((d - r0) / dr)
                            row_num = min(lig_elements.index(ele_lig), 22)
                            col_num = 99
                            local_number_array[i][row_num, col_num] += 1
        local_number_arrays_list.append(local_number_array)
    return lig_name_list, local_number_arrays_list



def drugresidue_score(rec_file, lig_file, potential_array, dr=0.2, r0=2.0, rc=6.0):
    lig_name_list, atom_number_arrays_list = get_contact_num_array(rec_file, lig_file, dr, r0, rc)
    
    score_list = []
    for n, lig_name in enumerate(lig_name_list):
        atom_number_array = atom_number_arrays_list[n]
        atom_energy_array = atom_number_array * potential_array
        potential = np.sum(atom_energy_array)
        score_list.append(potential)
    return lig_name_list, np.array(score_list, dtype=float)