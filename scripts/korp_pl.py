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
