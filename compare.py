from rdkit import Chem
import sys

sys.path.append("/home/leejinwon")
from acerxn import chem, process
import numpy as np
import compute_chg_and_bo_gurobi
import time


if __name__ == "__main__":
    smiles = sys.argv[1]
    # try:
    #    resolve = int(sys.argv[2])
    # except:
    #    resolve = 1
    # resolve = bool(resolve)
    molecule = chem.Molecule(smiles)
    (_, _, _, _, _, _, bond_mapping, _, _) = compute_chg_and_bo_gurobi.get_lists(
        molecule
    )
    print("bond mapping\n", bond_mapping)
    # print(resolve)
    # try:
    rd_original = molecule.get_rd_mol()
    chg_mol = molecule.get_chg()
    print("mol chg", chg_mol)
    molecule.adj_matrix = molecule.get_adj_matrix()
    molecule.atom_feature["chg"] = None
    molecule.bo_matrix = None
    t = time.time()
    # (
    #    chg_list,
    #    bo_matrix,
    #    chg_list2,
    #    bo_matrix2,
    # ) = compute_chg_and_bo_gurobi.compute_chg_and_bo_debug(
    #    molecule, chg_mol, resolve=resolve
    # )
    chg_list, bo_matrix = compute_chg_and_bo_gurobi.compute_chg_and_bo(
        molecule, chg_mol, HalogenConstraints=True, printOptLog=True, mode="fc"
    )
    actobo_time = time.time() - t
    print("hihi", list(zip(molecule.get_element_list(), chg_list)))
    molecule.atom_feature["chg"] = chg_list
    molecule.bo_matrix = bo_matrix
    rd_actobo = molecule.get_rd_mol()
    try:
        Chem.SanitizeMol(rd_actobo)
        print(
            "atom num same?",
            len(molecule.atom_feature["chg"])
            == len([atom for atom in rd_actobo.GetAtoms()]),
        )
        # rd_actobo = Chem.RemoveHs(rd_actobo)
        smi_actobo = Chem.MolToSmiles(rd_actobo)
        print(smi_actobo)
        print(smiles)
        print(smiles == smi_actobo)
        # print(Chem.MolToInchi(rd_original) == Chem.MolToInchi(rd_actobo))
        print(actobo_time)
    except:
        print("Sanitization Failed")

    # molecule.atom_feature["chg"] = chg_list2
    # molecule.bo_matrix = bo_matrix2
    # rd_actobo = molecule.get_rd_mol()
    # try:
    #    Chem.SanitizeMol(rd_actobo)
    #    smi_actobo = Chem.MolToSmiles(rd_actobo)
    #    print(smi_actobo)
    #    print(smiles == smi_actobo)
    #    # print(Chem.MolToInchi(rd_original) == Chem.MolToInchi(rd_actobo))
    #    print(actobo_time)
    # except:
    #    print("Sanitization Failed")
