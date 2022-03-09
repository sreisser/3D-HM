# 3D-HM
Calculate the 3D hydrophobic moment of a molecule 
based on the electrostatic potential on the surface as published 
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4052240/) .
Use either pqr as input or pdb, which will be converted 
into pqr.

## Installation

There are three different ways how to use 3D-HM. 
1) Amino acid sequence of a helical peptide as input.
In this case you need to install AmberTools to generate a pdb
model of the peptide. Please follow the instructions 
[here](https://ambermd.org/GetAmber.php)
2) pdb file as input. pdb2pqr is called with a given
force field to generate a pqr file (containing partial charges
and atomic radii). This is then used for surface and electrostatics
calculation.
3) pqr file as input. This gives you most control and flexibility.
You can also use special chemical groups or residues, or
nucleic acids in your structure.

### Create conda environment

```
conda create --name 3D-HM python=3.10
conda activate 3D-HM
```
Install requirements
```
pip install -r requirements.txt
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
DIR_3DHM = ''
```

Set APBS root directory:
```
DIR_APBS = ''
```

If you need to create helical pdb files with AmberTools,
set the AmberTools directory:
```
DIR_AMBERTOOLS = ''
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








