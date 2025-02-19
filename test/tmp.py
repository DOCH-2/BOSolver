from rdkit import Chem

m = Chem.MolFromXYZFile("./cobalt_ligand.xyz")
print(Chem.MolToMolBlock(m))
print(Chem.MolToXYZBlock(m))
