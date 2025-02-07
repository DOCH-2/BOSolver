import numpy as np
from rdkit import Chem


PT = Chem.GetPeriodicTable()


def get_lists(mol: Chem.Mol):
    # period, group
    period_list = np.array(
        [PT.GetGetPeriod(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
    )
    group_list = np.array([PT.GetGroup(atom.GetAtomicNum()) for atom in mol.GetAtoms()])

    # valence, atomic number
    z_list = np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    ve_list = np.zeros([2 if period_list[i] == 1 else 8 for i in range(len(z_list))])

    # adj, neighbor, bond, bond_mapping
    bonds = mol.GetBonds()
