#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  preprocess protein and ligand

"""

from utils.cov_util import *
from utils.vina_dock import *
import utils.opts as opts
import argparse
import pandas as pd
from rdkit.Chem.Scaffolds import MurckoScaffold
import os

def parse_args():
    """ parsing arguments """
    parser = argparse.ArgumentParser(description='preprocess.py',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    opts.preprocess_opts(parser)
    opt = parser.parse_args() 
    return opt
    
def main():
    opt = parse_args()
    # covert smiles to 3d rdkit mol
    mol = smi2_3d(opt.smiles)
    ligpdbqt = f"{opt.outdir}/temp_lig.pdbqt"
    # save mol to pdbqt file
    mol2pdbqt(mol, ligpdbqt)
    
    receptorpdbqt = opt.pdbqt
    outpdbqt = f"{opt.outdir}/docked_lig.pdbqt"
    # carried out docking
    center, box_size, num_modes = read_vina_config(opt.vinaconf)
    v = dock_vina(ligpdbqt, receptorpdbqt, center, box_size, outpdbqt)
    
    # # read docking result pdbqt and format to rdkit mol
    dockpdbqt2mol(outpdbqt, f"{opt.outdir}/docked_lig")
    
    # convert recepotr pdbqt to pdb
    pdbqttopdb(receptorpdbqt, f"{opt.outdir}/receptor.pdb")

if __name__ == "__main__":
    main()
