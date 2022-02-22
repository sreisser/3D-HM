# 3D-HM
Calculate the 3D hydrophobic moment of a molecule 
based on the electrostatic potential on the surface as published 
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4052240/) .
Use either pqr as input or pdb, which will be converted 
into pqr.

## Installation
### Use conda environment
```
conda create --name 3D-HM python=3.10
conda activate 3D-HM
```
Install pdb2pqr
```
pip install pdb2pqr
```


### APBS
NanoShaper is used for surface calculation, APBS for
electrostatics calculation. Both are included in the latest 
release of APBS. Download and precompiled binaries here:
https://github.com/Electrostatics/apbs/releases and unzip.

### Set config
Edit config.py:
Set 3D-HM root directory:
```
root_dir = ''
```

Set APBS root directory:
```
apbs_dir = ''
```

### Run program 
```
# basic with pdb as input, generates pqr first
python src/3D-HM.py peptide.pdb
# use directly pqr input
python src/3D-HM.py peptide.pqr
# more parameters to choose
python src/3D-HM.py -h
```








