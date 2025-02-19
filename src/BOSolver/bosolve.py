import pathlib
import argparse
import textwrap

from rdkit import Chem
import numpy as np

from compute_chg_and_bo import compute_chg_and_bo
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


class LineWidthFormatter(argparse.HelpFormatter):
    width = 90

    def __init__(self, prog):
        super().__init__(prog=prog, width=LineWidthFormatter.width)


def bosolve_argparser():
    parser = argparse.ArgumentParser(
        prog="bosolve",
        description="""Bond Perception (Bond Order and Formal Charge) Program 
        based on Integer Linear Programming""",
        formatter_class=LineWidthFormatter,
    )
    parser.add_argument(
        "xyz",
        type=str,
        help="""
        xyz file (path to text file containing coordinates in xyz format,
        should end with .xyz) or xyz block (text in xyz format) of a system""",
    )
    parser.add_argument("c", type=int, default=0, help="total charge of a system")

    output_group = parser.add_argument_group("output options")
    output_group.add_argument(
        "-o",
        "--output",
        type=str,
        default="smi",
        help="output format",
        choices=["smi", "mol"],
        required=True,
    )
    output_group.add_argument(
        "--silent",
        action="store_true",
        help="do not print the result to the standard output",
    )
    output_group.add_argument(
        "--saveto",
        type=str,
        help="save the result to a file (path to the file)",
    )

    bosolver_group = parser.add_argument_group("BOSolver options")

    bosolver_group.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="verbose mode. print BO solving progress and dump optimization log file",
    )
    bosolver_group.add_argument(
        "--noResolve", action="store_true", help="do not proceed resolve_chg step"
    )
    bosolver_group.add_argument(
        "--noCleanUp",
        dest="noCleanUpHeuristics",
        action="store_true",
        help="do not apply clean-up heuristics",
    )
    bosolver_group.add_argument(
        "-fc",
        "--fcmode",
        action="store_true",
        help="set formal charge separation minimization as the primary objective in optimize_bo step",
    )
    bosolver_group.add_argument(
        "--metalCenters",
        nargs="*",
        dest="MetalCenters",
        type=int,
        help="atom indices for metal centers",
        default=[],
    )
    return parser


def main(args):
    import sys

    if args.silent:
        sys.stdout = open("/dev/null", "w")

    xyzfileORblock = args.xyz
    if xyzfileORblock.endswith(".xyz"):
        mol = Chem.MolFromXYZFile(xyzfileORblock)
    else:
        mol = Chem.MolFromXYZBlock(xyzfileORblock)
    chg = args.c

    bosolver_args = {
        k: v for k, v in vars(args).items() if k not in ["xyz", "c", "output"]
    }

    mol = perceiveConn(mol)
    mol = assignBO(mol, chg, **bosolver_args)

    out_format = args.output
    if out_format == "smi":
        outtext = Chem.MolToSmiles(mol, allHsExplicit=True)
    elif out_format == "mol":
        outtext = Chem.MolToMolBlock(mol)
    else:
        raise ValueError(f"Invalid output format: {out_format}")

    print(outtext)

    if args.saveto is not None:
        if pathlib.Path(args.saveto).exists():
            print(f"File {args.saveto} already exists. Overwrite.")
        with open(pathlib.Path(args.saveto), "w") as f:
            f.write(outtext)

    if args.silent:
        sys.stdout = sys.__stdout__

    return outtext


def main_cli():
    main(bosolve_argparser().parse_args())


if __name__ == "__main__":
    args = bosolve_argparser().parse_args()
    main(args)
