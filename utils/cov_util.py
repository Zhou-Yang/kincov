
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import numpy as np
from pymol import cmd
import os
import pandas as pd

def add_si(smi):
	m = Chem.MolFromSmiles(smi)
#	react = ['[C;!$(C[#7]);H2,H3:1]>> [C:1][Si]', '[c;H:1]>> [C:1][Si]']
#	react = ['[C;!$(C[#7]);H2,H3:1]>> [C:1][Si]', '[c;H:1]>> [C:1][Si]', '[n;!$(NC=O);!H0:1]>>[n:1][Si]', '[N;!$(NC=O);!H0:1]>>[N:1][Si]']
	react = ['[C;!$(C[#7]);H2,H3:1]>> [C:1][Si]', '[n;!$(NC=O);!H0:1]>>[n:1][Si]', '[N;!$(NC=O);!H0:1]>>[N:1][Si]', '[C;$(CO);H3:1]>>[C:1][Si]']
	lm = []
	for r in react:
		rxn = AllChem.ReactionFromSmarts(r)
		ps = rxn.RunReactants((m,))
		for i in range(len(ps)):
			p0 = ps[i][0]
			Chem.SanitizeMol(p0)
			lm.append(Chem.MolToSmiles(p0))
	return lm


def alignsmi2pdb(smi, pdb, outpdb):
	m = Chem.MolFromSmiles(smi)
	ref = Chem.rdmolfiles.MolFromPDBFile(pdb)
	m2 = Chem.AddHs(m)
	confids = AllChem.EmbedMultipleConfs(m2, numConfs=50)  # 50 conformers
	
	mcs = rdFMCS.FindMCS([ref, m])
	patt = Chem.MolFromSmarts(mcs.smartsString)
	refMatch = ref.GetSubstructMatch(patt)
	mv = m2.GetSubstructMatch(patt)
	rmsds = []
	for i in range(len(confids)):
		rmsds.append(AllChem.AlignMol(m2, ref, i, atomMap=list(zip(mv, refMatch))))
	miniRMSD = AllChem.AlignMol(m2, ref, int(np.argmin(rmsds)), atomMap=list(zip(mv, refMatch)))
	Chem.MolToPDBFile(m2, outpdb, int(np.argmin(rmsds)))
	return miniRMSD


def si_distance(lm, protein_pdbqt, inputpdb_fname, cys_s_name, outname):
	smilist = []
	rmsdlist = []
	distlist = []
	namelist = []
	cmd.load(protein_pdbqt)
	for i in range(len(lm)):
		outpdb = f"{ outname }_{ i }.pdb"
		# align
		rmsd = alignsmi2pdb(lm[i], inputpdb_fname, outpdb)
		cmd.load(outpdb)
		
		smilist.append(lm[i])
		rmsdlist.append(rmsd)
		distlist.append(cmd.distance(f"name SI1 and { outname }_{ i }", cys_s_name))
		namelist.append(f"{ outname }_{ i }")
		os.remove(outpdb)
	return pd.DataFrame({'NAME': namelist, 'SMILES': smilist, 'RMSD': rmsdlist, 'DISTANCE': distlist})


def SitoWarhead(smi, warhead_smi):
	m = Chem.MolFromSmiles(smi)
#	react = ['[C:1][Si] >> [C:1]%s' % warhead_smi, '[c:1][Si] >> [c:1]%s' % warhead_smi]
#	react = ['[C:1][Si] >> [C:1]%s' % warhead_smi, '[c:1][Si] >> [c:1]%s' % warhead_smi, '[n:1][Si] >> [n:1]%s' % warhead_smi, '[N:1][Si] >> [N:1]%s' % warhead_smi ]
	react = ['[C:1][Si] >> [C:1]%s' % warhead_smi, '[n:1][Si] >> [n:1]%s' % warhead_smi, '[N:1][Si] >> [N:1]%s' % warhead_smi, '[O:1][C][Si] >> [O:1]%s' % warhead_smi ]
	lm = []
	for r in react:
		rxn = AllChem.ReactionFromSmarts(r)
		ps = rxn.RunReactants((m,))
		for i in range(len(ps)):
			p0 = ps[i][0]
			Chem.SanitizeMol(p0)
			lm.append(Chem.MolToSmiles(p0))
	return lm[0]


def sidf2wardf(si_df, frag_path):
    df = si_df
    distance_list = df['DISTANCE'].tolist()
    react = []
    w_d = []
    with open(frag_path, "r") as f:
        for line in f.readlines():
            react.append(line.strip('\n').split(',')[0])
            w_d.append(line.strip('\n').split(',')[1])
    smilist = []
    countID = 0
    idlist = []
    dist_list = []
    for i in range(len(df['SMILES'])):
    	for j in range(len(react)):
            if float(distance_list[i])-float(w_d[j]) > -2 and float(distance_list[i])-float(w_d[j]) < 2:
                smilist.append(SitoWarhead(df['SMILES'][i], react[j]))
                dist_list.append(float(distance_list[i])-float(w_d[j]))
                idlist.append(countID)
                countID += 1
            else:
                print("try %s,%s  distance not match" %(i,j))
    if countID == 0:
        print("no candidate generated")
    else:
        print("%s candidates generated" %(countID))
    return pd.DataFrame({'SMILES': smilist, 'NAME': idlist, 'DISTANCE': dist_list})


def sidf2wardf_simple(si_df, frag_path="./lib/fragment.smi"):
	df = si_df
	react = []
	with open(frag_path, "r") as f:
		for line in f.readlines():
			react.append(line.strip('\n'))
	smilist = []
	countID = 0
	idlist = []
	for i in range(len(df['SMILES'])):
		for warhead_smi in react:
			smilist.append(SitoWarhead(df['SMILES'][i], warhead_smi))
			idlist.append(countID)
			countID += 1
	war_dict = {'SMILES': smilist, 'NAME': idlist}
	df_war = pd.DataFrame(war_dict)
	return df_war


def pdbqttopdb(pdbqt, pdb):
	cmd.load(pdbqt)
	cmd.save(pdb)

