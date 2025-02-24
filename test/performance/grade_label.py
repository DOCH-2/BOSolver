"""test.py for label"""

import sys
import os.path as osp
from rdkit import Chem
from xyz2mol import read_xyz_file

from grade_noAns import check_All

fdir, chg, smi = sys.argv[1].strip('"').split()
#fdir, smi, chg = sys.argv[1].strip("'").strip('"').split()

atoms, _, _ = read_xyz_file(fdir)
natoms = len(atoms)

chg = int(chg)
name = osp.basename(fdir)

mol = Chem.MolFromSmiles(smi)
mol = Chem.AddHs(mol)
for atom in mol.GetAtoms():
    atom.SetNumRadicalElectrons(atom.GetTotalNumHs())
    atom.SetNumExplicitHs(0)
    atom.SetNoImplicit(True)
print(name, chg, *check_All(mol, chg), f"{0:.5f}", sep="\t")