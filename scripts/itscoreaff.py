import subprocess
import shutil
import uuid
import os


def itscoreaff_score(itscoreaff, pro_pdb, lig_mol2, tmp='./tmp'):
    os.makedirs(tmp, exist_ok=True)
    
    tmp_pro_mol2 = f'{tmp}/{uuid.uuid4().hex[:20]}.mol2'
    tmp_lig_mol2 = f'{tmp}/{uuid.uuid4().hex[:20]}.mol2'
    
    shutil.copy(lig_mol2, tmp_lig_mol2)
    cmd_obabel = f'obabel {pro_pdb} -O {tmp_pro_mol2}'

    try:
        obabel_result = subprocess.run(cmd_obabel, shell=True, capture_output=True, 
                                       text=True, check=True)
    
    except subprocess.CalledProcessError as e:
        print(f"ITScoreAff error while converting {pro_pdb} to mol2: {e}")
        return None

    cmd = f"{itscoreaff} {tmp_pro_mol2} {tmp_lig_mol2} | tail -n +2 | awk '{{print $NF}}'"
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        output_lines = result.stdout.strip().split('\n')
        os.remove(tmp_pro_mol2)
        os.remove(tmp_lig_mol2)
        return output_lines
        
    except subprocess.CalledProcessError as e:
        print(f"ITScoreAff Error: {e}")
        return None


if __name__ == '__main__':
    itscoreaff = '/data_test/home/lzzheng/wzh/softwares/ITScoreAff_v1.0/ITScoreAff'
    pro_pdb = '/data_test/home/lzzheng/wzh/softwares/statistical/integrated_model/examples/1a30_protein.pdb'
    lig_mol2 = '/data_test/home/lzzheng/wzh/softwares/statistical/integrated_model/examples/1a30_decoys.mol2'

    out = itscoreaff_score(itscoreaff, pro_pdb, lig_mol2)
    print(len(out))