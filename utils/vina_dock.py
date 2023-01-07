from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from vina import Vina
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
import pandas as pd

def read_vina_config(vina_config):
	a = pd.read_csv(vina_config, header=None, sep='=')
	center_x = float(a[a[0] == 'center_x '][1])
	center_y = float(a[a[0] == 'center_y '][1])
	center_z = float(a[a[0] == 'center_z '][1])
	size_x = float(a[a[0] == 'size_x '][1])
	size_y = float(a[a[0] == 'size_y '][1])
	size_z = float(a[a[0] == 'size_z '][1])
	center = [center_x, center_y, center_z]
	box_size = [size_x, size_y, size_z]
	num_modes = float(a[a[0] == 'num_modes '][1])
	return center, box_size, num_modes

def smi2_3d(smi):
	m = Chem.MolFromSmiles(smi)
	m2 = Chem.AddHs(m)
	AllChem.EmbedMolecule(m2)
	AllChem.MMFFOptimizeMolecule(m2)
	return m2


def mol2pdbqt(mol, out):
	preparator = MoleculePreparation()
	preparator.prepare(mol)
	preparator.write_pdbqt_file(out)
	return None


def dock_vina(lig, receptor, center, box_size, out):
	v = Vina(sf_name="vina", verbosity=1)
	v.set_receptor(receptor)
	v.set_ligand_from_file(lig)
	v.compute_vina_maps(center=center, box_size=box_size)
	v.dock(exhaustiveness=32, n_poses=20)
	v.write_poses(out, n_poses=20, overwrite=True)
	return None


def dockpdbqt2mol(pdbqt, dockoutpdb):
    pdbqt_mol = PDBQTMolecule.from_file(pdbqt, skip_typing=False)
    mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    for c in range(20):
        Chem.rdmolfiles.MolToPDBFile(mol[0],filename = f"{ dockoutpdb }{ c }.pdb", confId = c)
#    c = 0
#    for i in pdbqt_mol:
#	#   mol = i.export_rdkit_mol()
#        mol = RDKitMolCreate.from_pdbqt_mol(i)
#        Chem.MolToPDBFile(mol[0], f"{ dockoutpdb }{ c }.pdb")
##            print(Chem.MolToSmiles(mol[1]))
#        c += 1
#        print(c)
	#return pdbqt_mol[0].export_rdkit_mol()
    return RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)

