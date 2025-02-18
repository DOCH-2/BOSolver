from rdkit import Chem
from rdkit.Chem import AllChem

from compute_chg_and_bo_pulp import compute_chg_and_bo
from utils.coord2adj import BaberHodgkin, Simple, VdWRadius

print("main executed")

orimol = Chem.MolFromSmiles("C=C")
orimol = Chem.AddHs(orimol)
AllChem.EmbedMolecule(orimol)
print(AllChem.MMFFOptimizeMolecule(orimol))

xyz = Chem.MolToXYZBlock(orimol)
print(xyz)

mol = Chem.MolFromXYZBlock(xyz)
BPAlgo = Simple()
adj = BPAlgo(mol)
print(adj)
