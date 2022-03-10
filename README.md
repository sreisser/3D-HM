# 3D-HM
Calculate the 3D hydrophobic moment of a molecule 
based on the electrostatic potential on the surface as published 
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4052240/).

There are three different ways how to use 3D-HM. 
1) Amino acid sequence of a helical peptide as input.
In this case you need to install AmberTools to generate a pdb
model of the peptide (see below). 
2) pdb file as input. pdb2pqr is called with a given
force field to generate a pqr file (containing partial charges
and atomic radii). This is then used for surface and electrostatics
calculation.
3) pqr file as input. This gives you most control and flexibility.
You can also use special chemical groups or residues, or
nucleic acids in your structure.


## Installation

### AmberTools if starting from sequence
If you want to start from the sequence of a helical peptide, first
install AmberTools to generate a 3D pdb model of the peptide.
Please follow the instructions 
[here](https://ambermd.org/GetAmber.php).

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
release of APBS. Download precompiled binaries here:
https://github.com/Electrostatics/apbs/releases and unzip.

### Set config
Edit src/config.py:

```
DIR_3DHM = ''
DIR_APBS = ''
DIR_AMBERTOOLS = ''
```


### Run program 
```
# basic with pdb as input, generates pqr first
python src/3D-HM.py peptide.pdb
# for helical peptides: submit sequence only
# (compare examples/)
python src/3D-HM.py test.seq
# use directly pqr input
python src/3D-HM.py peptide.pqr
# more parameters to choose
python src/3D-HM.py -h
```

### References 
If you use 3D-HM for publications, please make sure to cite:

```
/* 3D-HM */
@article{reisser2014_3DHM,
    title = {3D Hydrophobic Moment Vectors as a Tool to Characterize the Surface Polarity of Amphiphilic Peptides},
    author = {Reisser, Sabine and Strandberg, Erik and Steinbrecher, Thomas and Ulrich, Anne S.},
    journal = {Biophysical Journal},
    volume = {106},
    number = {11},
    pages = {2385-2394},
    year = {2014},  
}

/* APBS and pdb2pqr */
@article{jurrus2018apbs,
    title = {Improvements to the APBS biomolecular solvation software suite: Improvements to the APBS Software Suite}, 
    author = {Jurrus, Elizabeth and Engel, Dave and Star, Keith and Monson, Kyle and Brandi, Juan and Felberg, Lisa E. and Brookes, David H. and Wilson, Leighton and Chen, Jiahui and Liles, Karina and Chun, Minju and Li, Peter and Gohara, David W. and Dolinsky, Todd and Konecny, Robert and Koes, David R. and Nielsen, Jens Erik and Head-Gordon, Teresa and Geng, Weihua and Krasny, Robert and Wei, Guo-Wei and Holst, Michael J. and McCammon, J. Andrew and Baker, Nathan A.}, 
    journal = {Protein Science}, 
    volume = {27}, 
    number = {1}, 
    pages = {112-128},
    year = {2018},
    publisher = {Wiley Blackwell (John Wiley & Sons)}, 
}

/* Nanoshaper */
@article{decherchi2013nanoshaper,
    title = {A general and Robust Ray-Casting-Based Algorithm for Triangulating Surfaces at the Nanoscale},
    author = {Decherchi, Sergio and Rocchia, Walter},
    journal = {PLoS ONE},
    year = {2013},
    month = {04},
    volume = {8},
    number = {4},
    pages = {e59744},
    publisher = {Public Library of Science},
} 


```


Please support the continued development of APBS by registering your use 
[here](http://eepurl.com/by4eQr)






