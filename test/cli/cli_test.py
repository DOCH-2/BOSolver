import sys

sys.path.append("../../src/BOSolver/")

from bosolve import main, bosolve_argparser

from rdkit import Chem
from rdkit.Chem import AllChem

COMMAND_STR = "CCO.xyz 0 -o smi"
COMMAND_OPTION = COMMAND_STR.split()

parser = bosolve_argparser()
args = parser.parse_args(COMMAND_OPTION)
main(args)
