# KinCov 
KinCov allow for the automatic generation of candidates from an electrophilic warhead fragment (EWF) library for covalent kinase inhibitor design.
## Dependencies
* Python (>=3.8)
* RDKit
* Pymol
* Meeko (optionally, for docking)

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
Then use [preprocess.py] for the docking process. 

```bash
python preprocess.py --smiles ${SMILES} --pdbqt ${PDBQT} --vinaconf $VINACONF --outdir ${OUTPUT}
```
This will generate docking poses for the next step. 

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


