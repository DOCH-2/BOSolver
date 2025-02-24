import pathlib
import argparse

from rdkit import Chem

from BOSolver.bosolve import perceiveConn, assignBO


class LineWidthFormatter(argparse.HelpFormatter):
    width = 90
    CLILOGO = "\n".join(
        [
            r"...........................................",
            r"______.._____._____......._................",
            r"|.___.\|.._../..___|.....|.|...............",
            r"|.|_/./|.|.|.\.`--...___.|.|_..._____._.__.",
            r"|.___.\|.|.|.|`--..\/._.\|.\.\././._.\.'__|",
            r"|.|_/./\.\_/./\__/./.(_).|.|\.V./..__/.|...",
            r"\____/..\___/\____/.\___/|_|.\_/.\___|_|...",
            r"...........................................",
            r"........................................... Original Work by KyungHoon Lee",
        ]
    )

    def __init__(self, prog):
        super().__init__(prog=prog, width=LineWidthFormatter.width)

    def format_help(self):
        help_text = super().format_help()
        return LineWidthFormatter.CLILOGO + "\n" * 2 + help_text


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

    if args.saveto is not None:
        if pathlib.Path(args.saveto).exists():
            print(f"File {args.saveto} already exists. Overwrite.")
        with open(pathlib.Path(args.saveto), "w") as f:
            f.write(outtext)

    if args.silent:
        sys.stdout = sys.__stdout__

    return outtext


def main_cli():
    parser = bosolve_argparser()
    print(main(parser.parse_args()))


if __name__ == "__main__":
    main_cli()
