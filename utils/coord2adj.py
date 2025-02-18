"""Coordinate-to-Adjacency Conversion (Bond Perception) Algorithms"""

import numpy as np
from scipy import spatial

from .chem import PT, Chem


class BondPerception:
    def __init__(self):
        pass

    def __call__(self, *args, **kwargs):
        raise NotImplementedError


class VdWRadius(BondPerception):
    """Van der Waals Radius-based Bond Perception Algorithm"""

    def __init__(self, relTol=1.15, absTol=0.45):
        self.relTol = relTol
        self.absTol = absTol

    def __call__(self, mol: Chem.Mol):
        n = mol.GetNumAtoms()
        radii = [PT.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
        radii_mat = np.repeat(radii, n).reshape(n, n)
        criteria_mat = (radii_mat + radii_mat.T) * self.relTol + self.absTol

        coords = mol.GetConformer().GetPositions()
        dist_mat = spatial.distance_matrix(coords, coords)
        adj = np.where(dist_mat < criteria_mat, 1, 0)
        np.fill_diagonal(adj, 0)

        return adj


class BaberHodgkin(VdWRadius):
    """Bond Perception Algorithm by Baber and Hodgkin
    Reference: J. Chem. Inf. Comput. Sci. 1992, 32, 401-406"""

    def __init__(self):
        self.relTol = 1.0
        self.absTol = 0.45
        super().__init__(relTol=self.relTol, absTol=self.absTol)


class Simple(VdWRadius):
    """Simple Bond Perception Algorithm
    $ \le \le $
    """

    def __init__(self):
        self.relTol = 1.15
        self.absTol = 0
        super().__init__(relTol=self.relTol, absTol=self.absTol)

    def __call__(self, mol: Chem.Mol):
        adj = super().__call__(mol)
        adj[(adj < 0.8) & (adj > 0)] = 1
        return adj
