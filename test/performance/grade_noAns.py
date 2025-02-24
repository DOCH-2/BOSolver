import sys
import os.path as osp

sys.path.append("/home/leejinwon")

import pickle
from copy import deepcopy

import numpy as np
from rdkit import Chem

group_number = {1: 1, 6: 4, 7: 5, 8: 6, 9: 7}


def check_boSum(mol: Chem.Mol):
    # mol = deepcopy(mol)
    Chem.KekulizeIfPossible(mol)
    # Chem.Kekulize(mol)
    bo = Chem.GetAdjacencyMatrix(mol, useBO=True)
    nbond = int(bo.sum() / 2)
    return nbond


def check_chgSepa(mol: Chem.Mol):
    formal_chgs = np.array([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    chg_sepa = np.sum(np.abs(formal_chgs))
    return chg_sepa


def check_nAromatic(mol: Chem.Mol):
    Chem.SanitizeMol(
        mol, catchErrors=True
    )  # if there is any error, it will be caught by check_sanitize
    aro_atom = [atom.GetIsAromatic() for atom in mol.GetAtoms()]
    nAro = sum(aro_atom)
    return nAro


def check_nOctet(mol: Chem.Mol):
    nOctet = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            ve = group_number[1] - atom.GetFormalCharge() + atom.GetTotalValence()
            if ve == 2:
                nOctet += 1
        elif atom.GetAtomicNum() in [6, 7, 8, 9]:
            ve = (
                group_number[atom.GetAtomicNum()]
                - atom.GetFormalCharge()
                + atom.GetTotalValence()
            )
            if ve == 8:
                nOctet += 1
    return nOctet


def check_sanitize(mol: Chem.Mol):
    flag = Chem.SanitizeMol(mol, catchErrors=True)
    return int(flag == 0)


def check_chgConserve(mol: Chem.Mol, chg: int):
    formal_chgs = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    return int(sum(formal_chgs) == chg)


def check_noRadical(mol: Chem.Mol):
    flag = Chem.SanitizeMol(
        mol, catchErrors=True
    )  # if there is any error, it will be caught by check_sanitize
    if flag == 0:
        for atom in mol.GetAtoms():
            atom.SetNumRadicalElectrons(atom.GetTotalNumHs())
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(True)
    radicals = [(atom.GetNumRadicalElectrons() % 2) for atom in mol.GetAtoms()]
    return int(sum(radicals) == 0)


def check_All(mol: Chem.Mol, chg: int):
    if mol is None:
        return [-1] * 7
    return [
        check_chgConserve(mol, chg),
        check_noRadical(mol),
        check_sanitize(mol),
        check_boSum(mol),
        check_chgSepa(mol),
        check_nAromatic(mol),
        check_nOctet(mol),
    ]


if __name__ == "__main__":
    name, smi, t = sys.argv[1].split()

    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol, explicitOnly=True)
    # chg = -1
    chg = 0
    t = float(t)
    print(
        name,
        check_chgConserve(mol, chg),
        check_noRadical(mol),
        check_nOctet(mol),
        check_boSum(mol),
        check_chgSepa(mol),
        check_nAromatic(mol),
        check_sanitize(mol),
        f"{t:.5f}",
        sep="\t",
    )
