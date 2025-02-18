import pathlib

from rdkit import Chem
import numpy as np

from compute_chg_and_bo_pulp import compute_chg_and_bo
from utils.coord2adj import BondPerception, CovalentRadius


def perceiveConn(
    mol: Chem.Mol, algorithm: type[BondPerception] = CovalentRadius, **kwargs
):
    """Perceive connectivity information of a molecule based on the coordinates of atoms.

    returns rdkit.Chem.Mol object with connectivity information.
    All bond orders are set to rdkit.Chem.BondType.UNSPECIFIED.
    """
    algo = algorithm(**kwargs)
    return algo(mol)


def assignBO(mol: Chem.Mol, chg: int, **kwargs):
    """assigns Bond Orders and formal charges to a molecule

    returns rdkit.Chem.Mol object with bond orders and formal charges (re-)assigned.
    """
    assert mol.GetNumBonds() != 0, (
        "No connectivity information is given. Do perceiveConn first."
    )

    resolve = not kwargs.get("noResolve", False)
    cleanup = not kwargs.get("noCleanUpHeuristics", False)

    chg_list, bo_matrix = compute_chg_and_bo(mol, chg, resolve=resolve, cleanup=cleanup)

    assert chg_list is not None and bo_matrix is not None, (
        "BOSolver failed to assign bond orders and formal charges."
    )

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


def main():
    import argparse

    parser = argparse.ArgumentParser(
        prog="bosolve",
        description="Bond Perception (Bond Order and Formal Charge) Program based on Integer Linear Programming",
    )
    parser.add_argument(
        "xyz",
        type=str,
        help="xyz file (path to text file containing coordinates in xyz format) or xyz block (text in xyz format) of a system",
    )
    parser.add_argument("c", type=int, default=0, help="total charge of a system")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="smi",
        help="output format",
        choices=["smi", "mol"],
    )
    parser.add_argument(
        "-s",
        "--saveto",
        type=str,
        help="save the result to a file (path to the file)",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="verbose mode. dump optimization log file",
    )
    parser.add_argument(
        "--noResolve", action="store_true", help="do not proceed resolve_chg step"
    )
    parser.add_argument(
        "--noCleanUp",
        dest="noCleanUpHeuristics",
        action="store_true",
        help="do not apply clean-up heuristics",
    )

    args = parser.parse_args()
    print(args)

    xyzfileORblock = args.xyz
    if xyzfileORblock.endswith(".xyz"):
        mol = Chem.MolFromMolFile(xyzfileORblock, sanitize=False, removeHs=False)
    else:
        mol = Chem.MolFromXYZBlock(xyzfileORblock, sanitize=False, removeHs=False)
    chg = args.c

    bosolver_args = {
        k: v for k, v in vars(args).items() if k not in ["xyz", "c", "output"]
    }

    mol = perceiveConn(mol)
    mol = assignBO(mol, chg, **bosolver_args)

    out_format = args.output
    if out_format == "smi":
        outtext = Chem.MolToSmiles(mol)
    elif out_format == "mol":
        outtext = Chem.MolToMolBlock(mol)
    else:
        raise ValueError(f"Invalid output format: {out_format}")

    if args.saveto is not None:
        if pathlib.Path(args.saveto).exists():
            print(f"File {args.saveto} already exists. Overwrite.")
        with open(pathlib.Path(args.saveto), "w") as f:
            f.write(outtext)


if __name__ == "__main__":
    main()
