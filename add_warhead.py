#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
  add warhead to a ligand

"""


from utils.cov_util import *
from utils.vina_dock import *
import utils.opts as opts
import argparse


def parse_args():
    """ parsing arguments """
    parser = argparse.ArgumentParser(description='add_warhead.py',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    opts.addwarhead_opts(parser)
    opt = parser.parse_args()
    return opt

def main():
    opt = parse_args()
    # covert smiles to 3d rdkit mol
    mol = smi2_3d(opt.smiles)
    # add Si to input SMILES
    smiles_si = add_si(opt.smiles)
    if opt.pose is not None:
        complex_structure = 'docked_lig{opt.pose}.pdb'
    elif opt.complex is not None:
        complex_structure = opt.complex 
    # measure the distance of Si to Cys
    si_df = si_distance(
                smiles_si,
                opt.receptor,
                complex_structure,
                f'resid { opt.site } and name SG',
                'temp'
     )
    # change Si to warhead
    sidf2wardf(si_df, opt.warhead).to_csv(f'candidate.csv', index=None)


if __name__ == "__main__":
    main()
