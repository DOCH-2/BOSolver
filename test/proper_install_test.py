from rdkit import Chem
from rdkit.Chem import AllChem
import sys

sys.path.append("../src/BOSolver")

import bosolve

params = Chem.SmilesParserParams()
params.removeHs = False


def smi2xyz2smi_test(smi, chg, additional_args=[]):
    mol = Chem.MolFromSmiles(smi, params=params)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    xyz = Chem.MolToXYZBlock(mol)

    parser = bosolve.bosolve_argparser()
    args = parser.parse_args([str(xyz), str(chg), "-o", "smi"] + additional_args)
    output = bosolve.main(args)

    mol_restored = Chem.MolFromSmiles(output, params=params)

    smi_ori = Chem.MolToSmiles(mol)
    smi_res = Chem.MolToSmiles(mol_restored)

    print("smi_ori:", smi_ori)
    print("smi_res:", smi_res)

    return smi_ori == smi_res


def test_1():
    smi = "CCCCCC"
    chg = 0

    assert smi2xyz2smi_test(smi, chg)


def test_2():
    smi = "C1CCCCC1"
    chg = 0

    assert smi2xyz2smi_test(smi, chg)


def test_3():
    """Formal Charge First Mode Test"""
    smi = "[H]c1c([H])c(N([H])N=C(Sc2nnnn2-c2c([H])c([H])c(C([H])([H])[H])c([H])c2[H])C(=O)[C+]([H])[H])c([H])c([H])c1Cl"
    chg = 1

    odmode = smi2xyz2smi_test(smi, chg)
    fcmode = smi2xyz2smi_test(smi, chg, additional_args=["-fc"])
    assert (odmode is False) and (fcmode is True)


def test_4():
    smi = "C12=CC=CC3=C1C(CC(P(C4=CC(C=CC=C5CC=C6)=C5C6=C4)C7=CC(C=CC=C8CC=C9)=C8C9=C7)=C3)=CC=C2"
    chg = 0

    assert smi2xyz2smi_test(smi, chg)


def test_5():
    """Metal Center Test"""
    smi = "[H][Co]([C]=O)([C]=O)([C]=O)[C]=O"
    chg = 0
    metal_centers = []

    xyzfile = "./cobalt_ligand.xyz"
    parser = bosolve.bosolve_argparser()
    args = parser.parse_args(
        [xyzfile, str(chg), "-o", "smi"]
        + ["-fc", "--metalCenter", " ".join(metal_centers)]
    )
    output = bosolve.main(args)

    mol_restored = Chem.MolFromSmiles(output, params=params)

    smi_res = Chem.MolToSmiles(mol_restored)

    print("smi ori:", smi)
    print("smi res:", smi_res)

    assert smi == smi_res
