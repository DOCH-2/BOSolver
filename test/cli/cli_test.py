import BOSolver
import pathlib
from BOSolver.cli import main, bosolve_argparser


def test_cli():
    pathtoclitest = pathlib.Path(__file__)
    COMMAND_STR = f"{pathtoclitest.parent / 'CCO.xyz'} 0 -o smi"
    COMMAND_OPTION = COMMAND_STR.split()

    parser = bosolve_argparser()
    args = parser.parse_args(COMMAND_OPTION)
    text = main(args)
    assert text == "[H][O][C]([H])([H])[C]([H])([H])[H]"
