# KinCov 
KinCov allow for the automatic generation of candidates from an electrophilic warhead fragment (EWF) library for covalent kinase inhibitor design.
## Dependencies
* Python (>=3.8)
* RDKit
* Pymol
* Meeko (optionally, for docking)
* Amber (>=18)
* GROMACS (>=2020.6)

## Installation
Create a new conda environment:

```bash
conda env create -f environment.yml
conda activate kincov
```

## Use
### Simple Usage
To run kincov, just edit the configuration in [run.sh](https://github.com/Zhou-Yang/kincov/edit/main/run.sh) and run the script
```bash
bash run.sh
```

### Generate covalent candidates without a co-crystal structure
Prepare a PDBQT file for the receptor. This can be prepared by Meeko.
see https://github.com/forlilab/Meeko.
Then use [preprocess.py](https://github.com/Zhou-Yang/kincov/edit/main/preprocess.py) for the docking process. 

```bash
python preprocess.py --smiles ${SMILES} --pdbqt ${PDBQT} --vinaconf $VINACONF --outdir ${OUTPUT}
```
This will generate docking poses for the next step. 

### Reactivity analysis of the candidate cysteines based on $pK_a$ values
The $pK_a$ values of the kinase cysteines will be predicted by CpHMD simulations using the method as described in [[1]](#1). This calculation is time-consuming and is recommended to be performed on a cluster. 

```bash
python pka.py --pdb ${PDB} --site ${SITE_ID} --outdir ${OUTPUT}
sbatch pka.sbatch 
```

### Accessibility of cysteines with high reactivity 
The position of a cysteine and its distance to a reversible inhibitor will be sampled using replica exchange with solute scaling (REST2) [[2]](#2). This calculation is time-consuming and is recommended to be performed on a cluster. 

```bash
python ac.py --pdb ${PDB} --site ${SITE_ID} --outdir ${OUTPUT}
sbatch ac.sbatch
```

### Add EWFs and estimate proximity of the EWF to the cysteine
The EWFs will be added automaticlly to the input reversible inhibitor
$$d_{esti}=d_{site}-L_{EWF}$$
where $d_{site}$ is the distance to the cysteine, and $L_{EWF}$ is the length of the EWF. 
```bash
python add_warhead.py --smiles ${SMILES} --receptor ${OUTPUT}/${RECEPTOR_PDB} --complex ${OUTPUT}/${COMPLEX_PDB} --outcsv ${OUTPUT}.csv --warhead ${WARHEAD_LIB} --site ${SITE_ID}
```
### Covalent docking (OPTIONAL)
Covalent docking can be carried out directly for the candidate generated
```bash
python covdock.py --smiles ${SMILES} --pdbqt ${OUTPUT}/${PDBQT} --site ${SITE_ID} --ligand ${OUTPUT}/ligand.pdb
```

# References
<a id="1">[1]</a> John Mongan, David A. Case, and J. Andrew McCammon, Constant pH molecular dynamics in generalized Born implicit solvent, J. Comput. Chem., 2004, 25 (16), 2038-2048. https://doi.org/10.1002/jcc.20139

<a id="2">[2]</a> Lingle Wang, Richard A. Friesner, and B. J. Berne. Replica Exchange with Solute Scaling: A More Efficient Version of Replica Exchange with Solute Tempering (REST2). J. Phys. Chem. B, 2011, 115 (30), 9431–9438. https://doi.org/10.1021/jp204407d
