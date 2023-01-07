from rdkit import Chem
from rdkit.Chem import AllChem
import argparse
import subprocess
from pymol import cmd

def addSC(local_smiles):
	lm = []
	ps = AllChem.ReactionFromSmarts('[N:1]C(C=C)=O >> [N:1]C(C=CSC)=O').RunReactants(
		(Chem.MolFromSmiles(local_smiles),))
	for local_i in range(len(ps)):
		p0 = ps[local_i][0]
		Chem.SanitizeMol(p0)
		lm.append(Chem.MolToSmiles(p0))
	return lm[0]

def getSCnumber(pdbfile):
	m = Chem.MolFromPDBFile(pdbfile)
	S_patt = Chem.MolFromSmarts('[S$(*C=CC=O)]')
	C_patt = Chem.MolFromSmarts('[C$(*SC=CC=O)]')
	if m.HasSubstructMatch(S_patt):
		localSid = m.GetSubstructMatch(S_patt)
	if m.HasSubstructMatch(C_patt):
		localCid = m.GetSubstructMatch(C_patt)
	return localSid[0], localCid[0]


def smi2pdb(local_smiles, pdbfile):
	m2 = Chem.AddHs(Chem.MolFromSmiles(local_smiles))
	AllChem.EmbedMolecule(m2)
	Chem.MolToPDBFile(m2, pdbfile)

def parse_args():
    """ parsing arguments """
    parser = argparse.ArgumentParser(description='covdock.py',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    opts.covdock(parser)
    opt = parser.parse_args()
    return opt

def main():
    opt = parse_args()
    smiles = args.smiles
    smi2pdb(smiles, f'ligand.pdb')
    site = f'A:CYS{opt.site}'
    outscpdb = f'ligand_SC.pdb'
    smi2pdb(addSC(smiles), outscpdb)
    S_id, C_id = getSCnumber(outscpdb)
    # kinase pdbqt
    receptor = "receptor.pdbqt"
    # mgl tools path
    MGL = '/home/kincov/AutoDock/mgltools_x86_64Linux2_1.5.7'
    # utilities24
    utilities24 = f'{mgl}/bin/pythonsh {mgl}/MGLToolsPckgs/AutoDockTools/Utilities24'
    execstring = f"""
    conda init bash
    source ~/.bashrc
    conda activate py27
    export ADCOV=/home/kincov/AutoDock/adCovalentDockResidue/adcovalent
    
    python $ADCOV/prepareCovalent.py --ligand {outscpdb} --ligindices {C_id + 1},{S_id + 1} --receptor {receptor} --residue {site} --outputfile ligcovalent.pdb
    
    {utilities24}/prepare_receptor4.py -r ligcovalent.pdb
    {utilities24}/prepare_flexreceptor4.py -r {receptor} -s receptor:{site}
    {utilities24}/prepare_flexreceptor4.py -r ligcovalent.pdbqt -s ligcovalent:{site}
    
    {utilities24}/prepare_gpf4.py -r receptor_rigid.pdbqt -x ligcovalent_flex.pdbqt -l ligcovalent_flex.pdbqt -o receptor.gpf -y -I 20
    {utilities24}/prepare_dpf4.py -r receptor_rigid.pdbqt -x ligcovalent_flex.pdbqt -l ligcovalent_flex.pdbqt -o ligcovalent_receptor.dpf -p move='empty'
    
    echo "rmsdatoms all" > empty
    sed -i 's/unbound_model extended/unbound_energy 0.0/g' ligcovalent_receptor.dpf
    
    export ATD=/home/kincov/AutoDock/adCovalentDockResidue/3upo/
    $ATD/autogrid4 -p receptor.gpf -l receptor.glg
    $ATD/autodock4 -p ligcovalent_receptor.dpf -l ligcovalent_receptor.dlg
    
    grep '^DOCKED' ligcovalent_receptor.dlg | cut -c9- > my_docking.pdbqt
    cut -c-66 my_docking.pdbqt > my_docking.pdb
    """
    subprocess.run(execstring, shell=True, executable='/bin/bash', check=True)
    
    cmd.load('my_docking.pdb')
    cmd.split_states('my_docking')
    for i in range(10):
    	i += 1
    	cmd.save(f'my_docking_{str(i).rjust(4, "0")}.pdb', f'my_docking_{str(i).rjust(4, "0")}')
    

if __name__ == "__main__":
    main()

