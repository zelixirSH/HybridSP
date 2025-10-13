import subprocess


def korp_score(korp, pro_pdb, lig_mol2):
    cmd = f"{korp} --receptor {pro_pdb} --ligand {lig_mol2} --mol2 | tail -n +14 | awk '{{print $NF}}'"
    
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        
        output_lines = result.stdout.strip().split('\n')
        return output_lines
        
    except subprocess.CalledProcessError as e:
        print(f"KORP-PL Error: {e}")
        return None


if __name__ == '__main__':
    korp = '/data_test/home/lzzheng/wzh/softwares/KORP-PL'
    pro_pdb = '/data_test/home/lzzheng/wzh/softwares/statistical/integrated_model/examples/1a30_protein.pdb'
    lig_mol2 = '/data_test/home/lzzheng/wzh/softwares/statistical/integrated_model/examples/1a30_decoys.mol2'
    out = korp_score(korp, pro_pdb, lig_mol2)

    print(out)