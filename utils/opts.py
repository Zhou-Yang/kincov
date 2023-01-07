import argparse

def preprocess_opts(parser):
    """ preprocess options """
    group = parser.add_argument_group('Data')
    group.add_argument("--smiles", help="input smiles")
    group.add_argument("--pdbid", help="kinase pdb ID")
    group.add_argument("--pdbqt", help="input pdbqt")
    group.add_argument("--vinaconf", help="input vinaconf")
    group.add_argument("--outdir", help="output dir")

def addwarhead_opts(parser):
    """ add warhead options """
    group = parser.add_argument_group('Data')
    parser.add_argument("--smiles", help="input smiles")
    parser.add_argument("--pose", help="docked_lig[0-20].pdb, choose an ID")
    parser.add_argument("--complex", help="docked_lig0.pdb, choose input complex structure")
    parser.add_argument("--site", help="cys resid")
    parser.add_argument("--outcsv", help="output file name")
    parser.add_argument("--warhead", help="warhead lib")
    parser.add_argument("--receptor", help="receptor.pdb")

def covdock_opts(parser):
    group = parser.add_argument_group('Data')
    parser.add_argument("--smiles", help="input smiles")
    parser.add_argument("--site", help="reaction site, e,g 427 ")
    group.add_argument("--pdbqt", help="input pdbqt")
    group.add_argument("--ligand", help="input ligand")
