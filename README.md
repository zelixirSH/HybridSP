# HybridSP
HybridSP is a novel hybrid statistical potential that combines distance-dependent atom-atom, atom-residue, and orientation-dependent interactions into a unified scoring function. 
It demonstrates exceptional accuracy in protein-ligand docking and virtual screening, rivaling or even surpassing state-of-the-art deep learning models.

## Installation
The following packages should be installed:

    numpy
    openbabel

HybridSP consists of three statistical potentials: ITScoreAff (distance-dependent atom-atom potential), DrugResidue<sub>W</sub> (distance-dependent atom-residue potential), and KORP-PL (orientation-dependent atom-residue potential). 
Among these, ITScoreAff and KORP-PL need to be downloaded separately, and the binary files `ITScoreAff` and `KORP-PL` should be placed in the "potentials" directory. 

The download links are as follows:

    ITScoreAff: http://huanglab.phys.hust.edu.cn/ITScoreAff/
    KORP-PL: https://team.inria.fr/nano-d/software/korp-pl/

## Usage
The protein should be provided in `pdb` format, and the ligand in `mol2` format. 

HybridSP functions as a post-scoring model for protein-ligand complexes; therefore, the ligand must first be properly docked into the binding pocket using a docking program such as AutoDock Vina.

HybridSP accommodates a variety of statistical potential models, including: 

  HybridSP<sub>dk</sub>
  HybridSP<sub>scr</sub>
  HybridSP<sub>bl</sub>
  DrugResidue
  DrugResidue<sub>W</sub>
  DrugResiGrp
  DrugResiGrp<sub>W</sub>
  DrugScoreRe
  DrugScoreRe<sub>W</sub>
  DrugScoreGrp
  DrugScoreGrp<sub>W</sub>

The characteristics of each model are detailed in the corresponding paper.
