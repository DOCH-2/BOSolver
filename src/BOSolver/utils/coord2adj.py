"""Coordinate-to-Adjacency Conversion (Bond Perception) Algorithms"""

import numpy as np
from scipy import spatial

from BOSolver.utils.chem import PT, Chem

_ACERXN_Rcov_TABLE = dict(
    [
        (1, 0.31),
        (3, 1.28),
        (4, 0.96),
        (5, 0.84),
        (6, 0.76),
        (7, 0.71),
        (8, 0.66),
        (9, 0.57),
        (11, 1.66),
        (12, 1.41),
        (13, 1.21),
        (14, 1.11),
        (15, 1.07),
        (16, 1.05),
        (17, 1.02),
        (18, 0.76),
        (19, 2.03),
        (20, 1.76),
        (22, 1.60),
        (26, 1.42),  # average of lowspin and highspin
        (27, 1.38),  # average of lowspin and highspin
        (28, 1.24),
        (35, 1.20),
        (45, 1.42),
        (46, 1.39),
        (53, 1.39),
    ]
)
# Reference: Daltan Trans., 2008, 2832-2838.


class DEFAULT_Rcov_TABLE(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.update(_ACERXN_Rcov_TABLE)

    def __getitem__(self, key):
        if key in self:
            return super().__getitem__(key)
        else:
            return PT.GetRcovalent(key)


class BondPerception:
    def __init__(self):
        pass

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class CovalentRadius(BondPerception):
    """Class for Covalent Radius-based Bond Perception Algorithm

    Covalent radii of elements referring to ACERXN. (original: Daltan Trans., 2008, 2832-2838.)
    If the radius of an element is not found in the table, the covalent radius from rdkit.Chem.PeriodicTable is used instead.
    Covalent radii of elements for bond perception can be manually adjusted and added, by modifying *attr* R_TABLE (works similarly to Dict).
    """

    _default_r_table = DEFAULT_Rcov_TABLE()

    def __init__(self, relTol=1.1, absTol=0.0):
        self.relTol = relTol
        self.absTol = absTol
        self.R_TABLE = CovalentRadius._default_r_table

    def determine_adj(self, mol: Chem.Mol):
        n = mol.GetNumAtoms()
        radii = [self.R_TABLE[atom.GetAtomicNum()] for atom in mol.GetAtoms()]
        radii_mat = np.repeat(radii, n).reshape(n, n)
        criteria_mat = (radii_mat + radii_mat.T) * self.relTol + self.absTol

        coords = mol.GetConformer().GetPositions()
        dist_mat = spatial.distance_matrix(coords, coords)
        adj = np.where(dist_mat < criteria_mat, 1, 0)
        np.fill_diagonal(adj, 0)

        return adj, dist_mat

    def __call__(self, mol: Chem.Mol, **kwargs) -> Chem.Mol:
        adj, _ = self.determine_adj(mol)
        i, j = np.nonzero(adj > 0)
        adj_list = np.vstack((i[i < j], j[i < j])).T.tolist()  # only i < j

        mol = Chem.RWMol(mol)

        # remove original bonds, if any
        original_bonds = list(
            map(lambda x: (x.GetBeginAtomIdx(), x.GetEndAtomIdx()), mol.GetBonds())
        )
        for bond in original_bonds:
            mol.RemoveBond(bond[0], bond[1])

        # add newly perceived bonds
        for bond in adj_list:
            mol.AddBond(bond[0], bond[1], Chem.BondType.UNSPECIFIED)

        mol = Chem.Mol(mol)

        return mol


class BaberHodgkin(CovalentRadius):
    """Bond Perception Algorithm by Baber and Hodgkin
    Reference: J. Chem. Inf. Comput. Sci. 1992, 32, 401-406"""

    def __init__(self):
        self.relTol = 1.0
        self.absTol = 0.45
        super().__init__(relTol=self.relTol, absTol=self.absTol)


class Simple(CovalentRadius):
    """Simple Bond Perception Algorithm"""

    def __init__(self):
        self.relTol = 1.15
        self.absTol = 0
        super().__init__(relTol=self.relTol, absTol=self.absTol)

    def determine_adj(self, mol: Chem.Mol):
        adj, dist_mat = super().determine_adj(mol)
        adj[(dist_mat < 0.8) & (adj > 0)] = 1

        return adj, dist_mat
