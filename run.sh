#!/usr/bin/env bash

######## Configure ##########

SMILES='O=C(C1CNCC1)CC2=CC3=C(C4=CC=NC=C4)NN=C3C=C2'
VINACONF=./examples/vina.conf
PDB_ID=4qtc
PDBQT=./examples/4qtc.prot.pdbqt
COMPLEX_PDB=docked_lig0.pdb
RECEPTOR_PDB=receptor.pdb
OUTPUT=./result
LIB_PATH=./lib
WARHEAD_LIB=${LIB_PATH}/warhead.smi
SITE_ID=702 #Cys resid
############################


#### Dock ligand to receptor (OPTIONAL) ####
if [ -z ${SMILES+x} ];then
    echo "SMILES is unset"
else
    if [ ! -f ${PDBQT} ];then
        echo "PDBQT file not found"
    else
        echo ${WARHEAD_LIB}
        python preprocess.py --smiles ${SMILES} --pdbqt ${PDBQT} --vinaconf $VINACONF --outdir ${OUTPUT}
    fi
fi

#### add warhead based on the complex structure ####
if [ ! -f ${OUTPUT}/${COMPLEX_PDB} ];then
    echo "complex pdb not exist"
else
    echo "load complex: ${COMPLEX_PDB}"
    echo "load receptor: ${RECEPTOR_PDB}"
    python add_warhead.py --smiles ${SMILES} --receptor ${OUTPUT}/${RECEPTOR_PDB} --complex ${OUTPUT}/${COMPLEX_PDB} --outcsv ${OUTPUT}.csv --warhead ${WARHEAD_LIB} --site ${SITE_ID} 
fi
#### covalent docking (OPTIONAL) ####
if [ ! -f $ADCOV/prepareCovalent.py ];then
    echo "prepareCovalent.py not found, please install"
else
    python covdock.py --smiles ${SMILES} --pdbqt ${OUTPUT}/${PDBQT} --site ${SITE_ID} --ligand ${OUTPUT}/ligand.pdb
fi
