from rdkit import Chem
import numpy as np

from BOSolver.compute_chg_and_bo import compute_chg_and_bo
from BOSolver.utils.coord2adj import BondPerception, CovalentRadius


def perceiveConn(
    mol: Chem.Mol, algorithm: type[BondPerception] = CovalentRadius(), **kwargs
) -> Chem.Mol:
    """Perceive connectivity information of a molecule based on the coordinates of atoms.

    returns rdkit.Chem.Mol object with connectivity information.
    All bond orders are set to rdkit.Chem.BondType.UNSPECIFIED.
    """
    algo = algorithm
    return algo(mol, **kwargs)


def assignBO(mol: Chem.Mol, chg: int, **kwargs) -> Chem.Mol:
    """assigns Bond Orders and formal charges to a molecule

    returns rdkit.Chem.Mol object with bond orders and formal charges (re-)assigned.
    """
    assert mol.GetNumBonds() != 0, (
        "No connectivity information is given. Do perceiveConn first."
    )

    resolve = not kwargs.get("noResolve", False)
    cleanup = not kwargs.get("noCleanUpHeuristics", False)
    verbose = kwargs.get("verbose", False)

    compute_chg_and_bo_kwargs = {
        k: v
        for k, v in kwargs.items()
        if k not in ["noResolve", "noCleanUpHeuristics", "verbose"]
    }
    compute_chg_and_bo_kwargs["printOptLog"] = verbose

    chg_list, bo_matrix = compute_chg_and_bo(
        mol, chg, resolve=resolve, cleanup=cleanup, **compute_chg_and_bo_kwargs
    )

    if chg_list is None and bo_matrix is None:
        raise RuntimeError("BOSolver failed to assign bond orders and formal charges.")

    assert chg_list is not None and bo_matrix is not None

    # modify bond orders and formal charges
    mol = Chem.RWMol(mol)

    # set Bond Order
    i, j = np.nonzero(bo_matrix > 0)
    bond_list = np.vstack((i[i < j], j[i < j])).T.tolist()
    for bond in bond_list:
        mol.RemoveBond(bond[0], bond[1])
        mol.AddBond(bond[0], bond[1], Chem.BondType.values[bo_matrix[bond[0], bond[1]]])

    # set Formal Charge
    chg_list = chg_list.tolist()
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(chg_list[atom.GetIdx()])

    return Chem.Mol(mol)
