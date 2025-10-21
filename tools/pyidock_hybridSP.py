#!/usr/bin/env python

import os, sys
import argparse
import subprocess as sp
import uuid
import shutil
import tempfile
import csv

# you may change this file path to your idock binary executable file
IDOCK_BINARY="idock"

ALLOWED_TERMS = ['receptor', 'ligand', 'out',
                 'center_x', 'center_y', 'center_z',
                 'size_x', 'size_y', 'size_z',
                 'threads', 'conformations', 'tasks'
                 ]

def argument():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", dest="config", default="vina.config", type=str,
                        help="Configuration file.")

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    return args


def generate_new_configs(config_inp_fpath, config_out_fpath, output_dpath):
    configs = {}
    new_configs = {}
    with open(config_inp_fpath) as lines:
        for l in lines:
            if len(l.split("=")) == 2:
                key = l.split("=")[0].strip()
                if key in ALLOWED_TERMS:
                    configs[key] = l.split("=")[1].strip("\n").strip()
                elif key == "cpu":
                    configs['threads'] = l.split("=")[1].strip("\n").strip()
                elif key == "exhaustiveness":
                    configs['tasks'] = l.split("=")[1].strip("\n").strip()
                elif key == "num_modes":
                    configs['conformations'] = l.split("=")[1].strip("\n").strip()

    print(configs)

    # now change output
    # idock output directory, not file path
    real_output_fpath = configs['out']
    configs['out'] = output_dpath

    with open(config_out_fpath, 'w') as tof:
        for key in configs.keys():
            tof.write("{} = {} \n".format(key, configs[key]))
    tof.close()

    return configs, real_output_fpath

def convert_pdbqt_to_mol2(pdbqt_file, mol2_file):
    """Convert PDBQT file to MOL2 format using Open Babel"""
    cmd = f"obabel -ipdbqt {pdbqt_file} -omol2 -O {mol2_file}"
    print(f"Converting PDBQT to MOL2: {cmd}")
    job = sp.Popen(cmd, shell=True)
    job.communicate()
    
    if os.path.exists(mol2_file):
        print(f"Successfully converted {pdbqt_file} to {mol2_file}")
        return True
    else:
        print(f"Failed to convert {pdbqt_file} to {mol2_file}")
        return False

def convert_pdbqt_to_pdb(pdbqt_file, pdb_file):
    """Convert PDBQT file to PDB format using Open Babel"""
    cmd = f"obabel -ipdbqt {pdbqt_file} -opdb -O {pdb_file}"
    print(f"Converting PDBQT to PDB: {cmd}")
    job = sp.Popen(cmd, shell=True)
    job.communicate()
    
    if os.path.exists(pdb_file):
        print(f"Successfully converted {pdbqt_file} to {pdb_file}")
        return True
    else:
        print(f"Failed to convert {pdbqt_file} to {pdb_file}")
        return False

def calculate_hybrid_score_via_subprocess(protein_pdb, ligand_mol2, method="HybridSPdk", output_csv=None):
    """Calculate hybrid score by calling score.py via subprocess"""
    
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    score_script = os.path.join(script_dir, "../score.py")
    
    if not os.path.exists(score_script):
        print(f"Error: score.py not found at {score_script}")
        return [], []
    
    # Create temporary output file if not provided
    if output_csv is None:
        output_csv = tempfile.NamedTemporaryFile(suffix='.csv', delete=False).name
    
    # Build the command
    cmd = [
        "python", score_script,
        "-r", protein_pdb,
        "-l", ligand_mol2,
        "-m", method,
        "-o", output_csv
    ]
    
    print(f"Running hybrid score calculation: {' '.join(cmd)}")
    
    try:
        # Run the scoring script
        result = sp.run(cmd, capture_output=True, text=True, check=True)
        print("Score calculation completed successfully")
        print(f"STDOUT: {result.stdout}")
        if result.stderr:
            print(f"STDERR: {result.stderr}")
        
        # Parse the output CSV file
        lig_names = []
        scores = []
        with open(output_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                lig_names.append(row['Lig_name'])
                scores.append(float(row['score']))
        
        print(f"Parsed {len(scores)} scores from {output_csv}")
        
        # Clean up temporary file if we created it
        if output_csv.startswith('/tmp'):
            os.remove(output_csv)
        
        return lig_names, scores
        
    except sp.CalledProcessError as e:
        print(f"Error running score.py: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        
        # Clean up temporary file if we created it
        if output_csv.startswith('/tmp') and os.path.exists(output_csv):
            os.remove(output_csv)
        
        return [], []
    
    except Exception as e:
        print(f"Unexpected error during score calculation: {e}")
        
        # Clean up temporary file if we created it
        if output_csv.startswith('/tmp') and os.path.exists(output_csv):
            os.remove(output_csv)
        
        return [], []

def add_scores_to_pdbqt(pdbqt_file, scores, method="HybridSPdk"):
    """Add hybrid scores to PDBQT file as REMARK lines"""
    temp_file = pdbqt_file + ".tmp"
    
    with open(pdbqt_file, 'r') as infile, open(temp_file, 'w') as outfile:
        model_count = -1
        current_model_lines = []
        
        for line in infile:
            current_model_lines.append(line)
            
            if line.startswith("MODEL"):
                model_count += 1
                
            elif line.startswith("ENDMDL"):
                # Write the current model
                for model_line in current_model_lines:
                    outfile.write(model_line)
                
                # Add the hybrid score remark if we have a score for this model
                if model_count < len(scores):
                    score_remark = f"REMARK HybridSP {method} Score: {scores[model_count]:.5f}\n"
                    outfile.write(score_remark)
                
                # Reset for next model
                current_model_lines = []
                
        # Handle any remaining lines
        for model_line in current_model_lines:
            outfile.write(model_line)
    
    # Replace original file with modified file
    shutil.move(temp_file, pdbqt_file)
    print(f"Added hybrid scores to {pdbqt_file}")

def run_idock(config_fpath):
    """Run idock and the post-scoring with HybridSP """
    # create a temporary path now
    tmp_token = str(uuid.uuid4().hex)[:8]
    temp_output_dpath = f"/tmp/idock_{tmp_token}"
    os.makedirs(temp_output_dpath, exist_ok=True)

    temp_config_fpath = os.path.join(temp_output_dpath, "idock.config")
    # create a new config file
    configs, real_output_fpath = generate_new_configs(config_fpath, temp_config_fpath, temp_output_dpath)

    # run the idock, a new file will be generated
    idock_temp_output_fpath = os.path.join(temp_output_dpath, os.path.basename(configs['ligand']))
    try:
        cmd = f"{IDOCK_BINARY} --config {temp_config_fpath}"
        print(f"Running cmd: {cmd}")
        job = sp.Popen(cmd, shell=True)
        job.communicate()
    except Exception as e:
        print(f"Running idock failed: {e}")
        return

    # check output file and cp it
    src = idock_temp_output_fpath
    dst = real_output_fpath
    
    if os.path.exists(src):
        # Step 1: Convert docking results to MOL2 format for scoring
        mol2_output = os.path.join(temp_output_dpath, "docking_results.mol2")
        if convert_pdbqt_to_mol2(src, mol2_output):
            # Step 2: Prepare protein file
            protein_file = configs['receptor']
            
            # Ensure protein is in PDB format
            if protein_file.endswith('.pdbqt'):
                protein_pdb = os.path.join(temp_output_dpath, "protein.pdb")
                if convert_pdbqt_to_pdb(protein_file, protein_pdb):
                    protein_file = protein_pdb
                else:
                    print(f"Warning: Could not convert protein file {protein_file} to PDB format")
                    # Try to use original file anyway
                    protein_pdb = protein_file.replace('.pdbqt', '.pdb')
                    if os.path.exists(protein_pdb):
                        protein_file = protein_pdb
                    else:
                        print(f"Error: No suitable protein PDB file found")
                        return
            
            # Step 3: Calculate hybrid scores using subprocess
            print(f"Calculating hybrid scores for {mol2_output} using protein {protein_file}")
            lig_names, hybrid_scores = calculate_hybrid_score_via_subprocess(protein_file, mol2_output, "HybridSPdk")
            
            if hybrid_scores:
                print(f"Calculated {len(hybrid_scores)} hybrid scores")
                
                # Step 4: Add scores to the original PDBQT file
                add_scores_to_pdbqt(src, hybrid_scores, "HybridSPdk")
            else:
                print("No hybrid scores calculated")
        
        # Step 5: Copy the modified file to final destination
        tofile = open(dst, 'w')
        with open(src, 'r') as lines:
            for l in lines:
                if "REMARK 921   NORMALIZED" in l:
                    score = float(l.split()[-2])
                    tofile.write("REMARK VINA RESULT:  {:.3f}   0.000  0.000 \n".format(score))
                tofile.write(l)
        tofile.close()
        
        print(f"Successfully generated output with hybrid scores: {dst}")
    else:
        print(f"Idock result {idock_temp_output_fpath} not found, exit now ...")

    try:
        shutil.rmtree(temp_output_dpath)
    except:
        print(f"removing directory {temp_output_dpath} failed")

if __name__ == "__main__":
    args = argument()
    run_idock(args.config)