from rdkit import Chem
from acerxn import chem
from acerxn import process
import numpy as np
import compute_chg_and_bo_gurobi
import time
import sys


if __name__ == "__main__":
    smiles = sys.argv[1]
    try:
        resolve = int(sys.argv[2])
    except:
        resolve = 1
    resolve = bool(resolve)
    molecule = chem.Molecule(smiles)
    # print(resolve)
    # try:
    chg_mol = molecule.get_chg()
    print("mol chg", chg_mol)
    molecule.adj_matrix = molecule.get_adj_matrix()
    molecule.atom_feature["chg"] = None
    molecule.bo_matrix = None
    t = time.time()
    (
        chg_list,
        bo_matrix,
        chg_list2,
        bo_matrix2,
    ) = compute_chg_and_bo_gurobi.compute_chg_and_bo(molecule, chg_mol, resolve=resolve)
    actobo_time = time.time() - t
    print("hihi", molecule.get_element_list())
    print("hihi", chg_list, len(chg_list))
    print("hihi2", chg_list2, len(chg_list2))
    molecule.atom_feature["chg"] = chg_list
    molecule.bo_matrix = bo_matrix
    rd_actobo = molecule.get_rd_mol()
    try:
        # Chem.SanitizeMol(rd_actobo)
        smi_actobo = Chem.MolToSmiles(rd_actobo)
        print(smi_actobo)
        print(smiles == smi_actobo)
        print(actobo_time)
    except:
        print("Sanitization Failed")

    molecule.atom_feature["chg"] = chg_list2
    molecule.bo_matrix = bo_matrix2
    rd_actobo = molecule.get_rd_mol()
    try:
        # Chem.SanitizeMol(rd_actobo)
        smi_actobo = Chem.MolToSmiles(rd_actobo)
        print(smi_actobo)
        print(smiles == smi_actobo)
        print(actobo_time)
    except:
        print("Sanitization Failed")
