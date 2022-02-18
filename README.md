# 3D-HM
Hydrophobic moment calculator

### Installation
### Requirements
#### NanoShaper
3D-HM uses a triangulated surface created by NanoShaper.
NanoShaper can be installed using conda-forge:
```
conda config --add channels conda-forge 
conda config --set channel_priority strict
conda install -c electrostatics nanoshaper
```

Source code is available at https://gitlab.iit.it/SDecherchi/nanoshaper

### APBS
For electrostatics calculation, get latest release here:
https://github.com/Electrostatics/apbs/releases

#### ChimeraX
Get ChimeraX from github: https://github.com/RBVI/ChimeraX.git
```
cd 3D-HM
git clone https://github.com/RBVI/ChimeraX.git
# link surface calculation files 
ln -s ChimeraX/src/bundles/surface/src/sasa.py src/
```






