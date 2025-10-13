import sys
import os
import argparse
import numpy as np


script_dir = os.path.dirname(os.path.abspath(__file__))
korp = os.path.join(f'{script_dir}/potentials/', 'KORP-PL')
itscoreaff = os.path.join(f'{script_dir}/potentials/', 'ITScoreAff')

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--receptor', type=str, default="examples/1a30_protein.pdb", 
                    help="Input protein pdb file.")
parser.add_argument('-l', '--ligand', type=str, default="examples/1a30_decoys.mol2", 
                    help="Input ligand mol2 file.")
parser.add_argument('-o', '--output', type=str, default="score.txt",
                    help="Output score file.")
parser.add_argument('-m', '--method', type=str, default="hybridSPdk",
                    choices=["HybridSPdk", "HybridSPscr", "HybridSPbl", 
                             "DrugResidue", "DrugResidue_W", "DrugResiGrp",
                             "DrugResiGrp_W", "DrugScoreRe", "DrugScoreRe_W", 
                             "DrugScoreGrp", "DrugScoreGrp_W"], 
                    help="Choose the statistical potential.")

args = parser.parse_args()

if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(0)

lig_mol2 = args.ligand
rec_pdb = args.receptor

if not lig_mol2.endswith('.mol2'):
    print(f'Error: The input ligand file {lig_mol2} should be converted into .mol2 format.')
    sys.exit(0)

if not rec_pdb.endswith('.pdb'):
    print(f'Error: The input protein file {rec_pdb} should be converted into .pdb format.')
    sys.exit(0)


if args.method == "HybridSPdk":
    from scripts.drugresidue import drugresidue_score
    from scripts.korp_pl import korp_score
    from scripts.itscoreaff import itscoreaff_score

    tmp_dir = './tmp'
    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugResidue_W.npy'))
    lig_name_list, score_drugresidue_w_array = drugresidue_score(rec_pdb, lig_mol2, potential_array)
    score_itscoreaff_array = np.array(itscoreaff_score(itscoreaff, rec_pdb, lig_mol2, tmp_dir), dtype=float)
    score_korp_array = np.array(korp_score(korp, rec_pdb, lig_mol2), dtype=float)
    os.rmdir(tmp_dir)

    a = 0.4
    b = 0.5
    c = 0.1
    score_arr = a*score_drugresidue_w_array + b*score_itscoreaff_array + c*score_korp_array


elif args.method == "HybridSPscr":
    from scripts.drugresidue import drugresidue_score
    from scripts.korp_pl import korp_score
    from scripts.itscoreaff import itscoreaff_score

    tmp_dir = './tmp'
    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugResidue_W.npy'))
    lig_name_list, score_drugresidue_w_array = drugresidue_score(rec_pdb, lig_mol2, potential_array)
    score_itscoreaff_array = np.array(itscoreaff_score(itscoreaff, rec_pdb, lig_mol2, tmp_dir), dtype=float)
    score_korp_array = np.array(korp_score(korp, rec_pdb, lig_mol2), dtype=float)
    os.rmdir(tmp_dir)
    
    a = 0.6
    b = 0.2
    c = 0.2
    score_arr = a*score_drugresidue_w_array + b*score_itscoreaff_array + c*score_korp_array


elif args.method == "HybridSPbl":
    from scripts.drugresidue import drugresidue_score
    from scripts.korp_pl import korp_score
    from scripts.itscoreaff import itscoreaff_score

    tmp_dir = './tmp'
    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugResidue_W.npy'))
    lig_name_list, score_drugresidue_w_array = drugresidue_score(rec_pdb, lig_mol2, potential_array)
    score_itscoreaff_array = np.array(itscoreaff_score(itscoreaff, rec_pdb, lig_mol2, tmp_dir), dtype=float)
    score_korp_array = np.array(korp_score(korp, rec_pdb, lig_mol2), dtype=float)
    os.rmdir(tmp_dir)

    a = 0.5
    b = 0.3
    c = 0.2
    score_arr = a*score_drugresidue_w_array + b*score_itscoreaff_array + c*score_korp_array


elif args.method == "DrugResidue":
    from scripts.drugresidue import drugresidue_score

    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugResidue.npy'))
    lig_name_list, score_drugresidue_array = drugresidue_score(rec_pdb, lig_mol2, potential_array)
    score_arr = score_drugresidue_array


elif args.method == "DrugResidue_W":
    from scripts.drugresidue import drugresidue_score

    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugResidue_W.npy'))
    lig_name_list, score_drugresidue_w_array = drugresidue_score(rec_pdb, lig_mol2, potential_array)
    score_arr = score_drugresidue_w_array


elif args.method == "DrugResiGrp":
    from scripts.drugresidue import drugresidue_score
    from scripts.group import group_score
    
    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugResidue.npy'))
    lig_name_list, score_drugresidue_array = drugresidue_score(rec_pdb, lig_mol2, potential_array)
    grp_potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_group.npy'))
    score_group_array = group_score(rec_pdb, lig_mol2, grp_potential_array)
    score_arr = 0.9*score_drugresidue_array + 0.1*score_group_array


elif args.method == "DrugResiGrp_W":
    from scripts.drugresidue import drugresidue_score
    from scripts.group import group_score

    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugResidue_W.npy'))
    lig_name_list, score_drugresidue_w_array = drugresidue_score(rec_pdb, lig_mol2, potential_array)
    grp_potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_group_W.npy'))
    score_group_w_array = group_score(rec_pdb, lig_mol2, grp_potential_array)
    score_arr = 0.9*score_drugresidue_w_array + 0.1*score_group_w_array


elif args.method == "DrugScoreRe":
    from scripts.drugscore import drugscore_score

    tmp_dir = './tmp'
    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugScoreRe.npy'))
    lig_name_list, score_drugscore_array = drugscore_score(rec_pdb, lig_mol2, potential_array, tmp_dir)
    score_arr = score_drugscore_array
    os.rmdir(tmp_dir)


elif args.method == "DrugScoreRe_W":
    from scripts.drugscore import drugscore_score

    tmp_dir = './tmp'
    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugScoreRe_W.npy'))
    lig_name_list, score_drugscore_w_array = drugscore_score(rec_pdb, lig_mol2, potential_array, tmp_dir)
    score_arr = score_drugscore_w_array
    os.rmdir(tmp_dir)


elif args.method == "DrugScoreGrp":
    from scripts.drugscore import drugscore_score
    from scripts.group import group_score

    tmp_dir = './tmp'
    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugScoreRe.npy'))
    lig_name_list, score_drugscore_array = drugscore_score(rec_pdb, lig_mol2, potential_array, tmp_dir)
    grp_potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_group.npy'))
    score_group_array = group_score(f"{rec_pdb.split('.')[0]}.pdb", lig_mol2, grp_potential_array)
    score_arr = 0.9*score_drugscore_array + 0.1*score_group_array
    os.rmdir(tmp_dir)


elif args.method == "DrugScoreGrp_W":
    from scripts.drugscore import drugscore_score
    from scripts.group import group_score

    tmp_dir = './tmp'
    potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_DrugScoreRe_W.npy'))
    lig_name_list, score_drugscore_w_array = drugscore_score(rec_pdb, lig_mol2, potential_array, tmp_dir)
    grp_potential_array = np.load(os.path.join(f'{script_dir}/potentials/', 'potential_group_W.npy'))
    score_group_w_array = group_score(f"{rec_pdb.split('.')[0]}.pdb", lig_mol2, grp_potential_array)
    score_arr = 0.9*score_drugscore_w_array + 0.1*score_group_w_array
    os.rmdir(tmp_dir)



for n in range(len(lig_name_list)):
    print(n, lig_name_list[n], score_arr[n])


with open(args.output, 'w') as f:
    f.write('Lig_ID,Lig_name,score\n')
    for n in range(len(lig_name_list)):
        f.write(f"{n+1},{lig_name_list[n]},{score_arr[n]:.5f}\n")
