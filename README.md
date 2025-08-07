# CTD Clustering
## Clustering of RNA Polymerase II C-Terminal Domain Models using MD simulations
This repository contains the scripts for the analysis of the simulations of concentrated 2CTD (CTD with 2-heptapeptide repeats) non-phosphorylated and phosphorylated systems. It also contains extracted representative frames with the lowest energy conformations according to the PCA of the cartesian coordinates of the peptides along the simulation trajectories.

*** Requirements:
```
Python versions 3+

MDAnalysis

tidynamics
```
*** Usage:
```
Python [options] filename.py
```
*** Example:

To calculate inter protein-protein contacts:
```
python contacts.py --traj_path=dcdfile --psf_path=psffile --sela="protein" --selb="protein" --out_path=outputfile
```
*** Citation:
```
Amith, W. D. Chen, V. T. Dutagaci B., Clustering of RNA Polymerase II C-Terminal Domain Models upon Phosphorylation, Journal of Physical Chemistry B 2024, 128, 10385−10396
