from scripts.util import lig_elements, resi_elements_dic, drugscore_elements
import numpy as np


def readmol2(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]

    molecules = []
    current_name = None
    current_sybyl = []
    current_xyz = []
    
    all_xyz = []
    inside_atom_section = False
    
    for i, line in enumerate(lines):
        if line.startswith('@<TRIPOS>MOLECULE'):
            # 如果在一个分子的过程中进入新分子，则保存当前分子数据
            if current_name is not None:
                current_sybyl = np.array(current_sybyl)
                current_xyz = np.array(current_xyz).astype(float)
                molecules.append((current_name, current_sybyl, current_xyz))
            
            # 开始新的分子
            current_name = lines[i + 1]  # 分子名在下一行
            current_sybyl = []
            current_xyz = []
            inside_atom_section = False
        
        elif line.startswith('@<TRIPOS>ATOM'):
            inside_atom_section = True  # 标记进入 ATOM 区域
        
        elif line.startswith('@<TRIPOS>BOND'):
            inside_atom_section = False  # 标记离开 ATOM 区域
        
        elif inside_atom_section and line:  # 在 ATOM 区域内解析原子信息
            atom_data = line.split()
            if len(atom_data) >= 6:  # 确保原子信息完整
                atom_name = atom_data[5]
                if atom_name in lig_elements:
                    current_sybyl.append(atom_name)
                    current_xyz.append(atom_data[2:5])
    
    if current_name is not None:
        current_sybyl = np.array(current_sybyl)
        current_xyz = np.array(current_xyz).astype(float)
        molecules.append((current_name, current_sybyl, current_xyz))
        all_xyz.extend(current_xyz)
    
    all_xyz = np.array(all_xyz).astype(float)
    min_xyz = np.min(all_xyz, axis=0) - 20.0
    max_xyz = np.max(all_xyz, axis=0) + 20.0
    return molecules, min_xyz, max_xyz


def readpdb(filename, min_xyz, max_xyz):
    min_x, min_y, min_z = min_xyz
    max_x, max_y, max_z = max_xyz

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
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if min_x < x < max_x and min_y < y < max_y and min_z < z < max_z:
                    if (residue in resi_elements_dic and atom in resi_elements_dic[residue]) or \
                       (residue in resi_elements_dic['OTH']):
                        data_atom_resi.append([atom, residue])
                        data_x.append(x)
                        data_y.append(y)
                        data_z.append(z)
                    
    data_atom_resi = np.array(data_atom_resi)
    data_xyz = list(zip(data_x, data_y, data_z))
    data_xyz = np.array(data_xyz)
    return data_atom_resi, data_xyz


def read_lig_mol2_drugscore(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]

    molecules = []
    current_name = None
    current_sybyl = []
    current_xyz = []
    
    all_xyz = []
    inside_atom_section = False
    
    for i, line in enumerate(lines):
        if line.startswith('@<TRIPOS>MOLECULE'):
            # 如果在一个分子的过程中进入新分子，则保存当前分子数据
            if current_name is not None:
                current_sybyl = np.array(current_sybyl)
                current_xyz = np.array(current_xyz).astype(float)
                molecules.append((current_name, current_sybyl, current_xyz))
            
            # 开始新的分子
            current_name = lines[i + 1]  # 分子名在下一行
            current_sybyl = []
            current_xyz = []
            inside_atom_section = False
        
        elif line.startswith('@<TRIPOS>ATOM'):
            inside_atom_section = True  # 标记进入 ATOM 区域
        
        elif line.startswith('@<TRIPOS>BOND'):
            inside_atom_section = False  # 标记离开 ATOM 区域
        
        elif inside_atom_section and line:  # 在 ATOM 区域内解析原子信息
            atom_data = line.split()
            if len(atom_data) >= 6:  # 确保原子信息完整
                atom_name = atom_data[5]
                if atom_name in drugscore_elements:
                    current_sybyl.append(atom_name)
                    current_xyz.append(atom_data[2:5])
    
    if current_name is not None:
        current_sybyl = np.array(current_sybyl)
        current_xyz = np.array(current_xyz).astype(float)
        molecules.append((current_name, current_sybyl, current_xyz))
        all_xyz.extend(current_xyz)
    
    all_xyz = np.array(all_xyz).astype(float)
    min_xyz = np.min(all_xyz, axis=0) - 20.0
    max_xyz = np.max(all_xyz, axis=0) + 20.0
    return molecules, min_xyz, max_xyz


def read_rec_mol2_drugscore(filename):
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
            residue = line.split()[-2][:3]
            if residue != 'HOH':
                data_sybyl.append(line.split()[5])
                data_xyz.append(line.split()[2:5])

    data_sybyl = np.array(data_sybyl)
    data_xyz = np.array(data_xyz).astype(float)
    return data_sybyl, data_xyz