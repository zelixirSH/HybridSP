# HybridSP
HybridSP is a novel hybrid statistical potential that combines distance-dependent atom-atom, atom-residue, and orientation-dependent interactions into a unified scoring function. It demonstrates exceptional accuracy in protein-ligand docking and virtual screening, rivaling or even surpassing state-of-the-art deep learning models.

## Installation
The following packages should be installed:

    numpy
    openbabel

HybridSP consists of three statistical potentials: ITScoreAff (distance-dependent atom-atom potential), DrugResidue$_\text{W}$ (distance-dependent atom-residue potential), and KORP-PL (orientation-dependent atom-residue potential). Among these, ITScoreAff and KORP-PL need to be downloaded separately, and the binary files `ITScoreAff` and `KORP-PL` should be placed in the "potentials" directory. The download links are as follows:
    ITScoreAff: http://huanglab.phys.hust.edu.cn/ITScoreAff/
    KORP-PL: https://team.inria.fr/nano-d/software/korp-pl/
