from rdkit import Chem
from rdkit.Chem import AllChem
import sys

sys.path.append("../src/BOSolver")

from utils import chem, coord2adj
import compute_chg_and_bo
import bosolve

mol = Chem.MolFromSmiles("CC")
mol = Chem.AddHs(mol)
# print(chem.get_lists(mol))

for bond in mol.GetBonds():
    print(
        bond.GetIdx(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType()
    )

print("-------------")

AllChem.EmbedMolecule(mol)
AllChem.MMFFOptimizeMolecule(mol)
mol = Chem.MolFromXYZBlock(Chem.MolToXYZBlock(mol))
mol = Chem.RWMol(mol)
mol.AddBond(0, 1, Chem.BondType.SINGLE)


BPAlgo = coord2adj.CovalentRadius()
mol = BPAlgo(mol)

for bond in mol.GetBonds():
    print(
        bond.GetIdx(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType()
    )

print("-------------")

chg_list, bo_matrix = compute_chg_and_bo.compute_chg_and_bo(mol, 0)

print(chg_list)
print(bo_matrix)

mol = Chem.MolFromSmiles("C1=CC=C2C=C3C=CC=CC3=CC2=C1")
mol = Chem.AddHs(mol)
from pprint import pprint

atoms_in_ring, bonds_in_ring, ring_neighbors_info = chem.get_ring_info(mol)
print(Chem.MolToSmiles(mol))
print()
pprint(atoms_in_ring)
print()
pprint(bonds_in_ring)
print()
pprint(ring_neighbors_info)

print()
print("-------------")

print(Chem.MolToSmiles(mol))
mol = bosolve.assignBO(mol, 0)
print(Chem.MolToSmiles(mol))
