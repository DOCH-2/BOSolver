from rdkit import Chem
from acerxn import chem
from acerxn import process
import numpy as np
import compute_chg_and_bo_gurobi
import time
import sys

if __name__ == '__main__':
    smiles = sys.argv[1]
    try:
        resolve = int(sys.argv[2])
    except:
        resolve = 1
    resolve = bool(resolve)
    molecule = chem.Molecule(smiles)
    print (resolve)
    #try:
    chg_mol = molecule.get_chg()
    molecule.adj_matrix = molecule.get_adj_matrix()
    molecule.atom_feature['chg'] = None
    molecule.bo_matrix = None
    t = time.time()
    chg_list, bo_matrix = compute_chg_and_bo_gurobi.compute_chg_and_bo(molecule, chg_mol, resolve=resolve)
    actobo_time = time.time()-t
    print(molecule.get_element_list())
    print(chg_list)
    molecule.atom_feature['chg'] = chg_list
    molecule.bo_matrix = bo_matrix
    rd_actobo = molecule.get_rd_mol()
    try:
        Chem.SanitizeMol(rd_actobo)
        smi_actobo = Chem.MolToSmiles(rd_actobo)
        print(smi_actobo)
        print(actobo_time)
    except:
        print('Sanitization Failed')
