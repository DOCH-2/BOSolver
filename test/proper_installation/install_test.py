import pathlib
from sys import path

from rdkit import Chem
from rdkit.Chem import AllChem

from BOSolver import bosolve

params = Chem.SmilesParserParams()
params.removeHs = False


def smi2xyz2smi_test(smi, chg, **kwargs):
    mol = Chem.MolFromSmiles(smi, params=params)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    xyz = Chem.MolToXYZBlock(mol)
    mol_restored = Chem.MolFromXYZBlock(xyz)
    mol_restored = bosolve.perceiveConn(mol_restored)
    mol_restored = bosolve.assignBO(mol_restored, chg, **kwargs)
    Chem.SanitizeMol(mol_restored)

    smi_ori = Chem.MolToSmiles(mol)
    smi_res = Chem.MolToSmiles(mol_restored)

    # print("smi_ori:", smi_ori) # for debugging
    # print("smi_res:", smi_res) # for debugging

    return smi_ori == smi_res


def xyz2smi_test(smi, xyzfile, chg, **kwargs):
    mol = Chem.MolFromXYZFile(xyzfile)
    mol = Chem.MolFromXYZFile(xyzfile)
    mol = bosolve.perceiveConn(mol)
    mol = bosolve.assignBO(mol, chg, **kwargs)
    Chem.SanitizeMol(mol)

    smi_res = Chem.MolToSmiles(mol)
    print("smi_res:", smi_res)  # for debugging

    return smi == smi_res


def test_1():
    smi = "CCCCCC"
    chg = 0

    assert smi2xyz2smi_test(smi, chg)


def test_2():
    smi = "C1CCCCC1"
    chg = 0

    assert smi2xyz2smi_test(smi, chg)


def test_3():
    """Formal Charge First Mode Test

    This example shows that occasionally people prefer minimizing the formal charge separation to making all atoms satisfy the octet rule.
    """
    smi = "[H]c1c([H])c(N([H])N=C(Sc2nnnn2-c2c([H])c([H])c(C([H])([H])[H])c([H])c2[H])C(=O)[C+]([H])[H])c([H])c([H])c1Cl"  # Carbon is not octet!
    chg = 1

    mol = Chem.MolFromSmiles(smi, params=params)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    xyz = Chem.MolToXYZBlock(mol)

    mol_od = bosolve.perceiveConn(Chem.MolFromXYZBlock(xyz))
    mol_od = bosolve.assignBO(mol_od, chg)
    Chem.SanitizeMol(mol_od)
    odmode = smi == Chem.MolToSmiles(mol_od)

    mol_fc = bosolve.perceiveConn(Chem.MolFromXYZBlock(xyz))
    mol_fc = bosolve.assignBO(mol_fc, chg, fcmode=True)
    Chem.SanitizeMol(mol_fc)
    fcmode = smi == Chem.MolToSmiles(mol_fc)

    assert (odmode is False) and (fcmode is True)


def test_4():
    smi = "C12=CC=CC3=C1C(CC(P(C4=CC(C=CC=C5CC=C6)=C5C6=C4)C7=CC(C=CC=C8CC=C9)=C8C9=C7)=C3)=CC=C2"
    chg = 0

    assert smi2xyz2smi_test(smi, chg)


def test_5():
    """Metal Center Test 1"""
    smi = "[H]C1=c2c(C([H])([H])C([H])([H])[H])c(C([H])([H])C([H])([H])[H])c3n2[Pt]24n5c1c(C([H])([H])C([H])([H])[H])c(C([H])([H])C([H])([H])[H])c5C([H])=c1c(C([H])([H])C([H])([H])[H])c(C([H])([H])C([H])([H])[H])c(n12)=C([H])c1c(C([H])([H])C([H])([H])[H])c(C([H])([H])C([H])([H])[H])c(n14)C=3[H]"
    chg = 0
    metal_centers = [0]

    xyzfile = (
        (pathlib.Path(__file__).parent / "./platinum_ligand.xyz").absolute().as_posix()
    )
    assert xyz2smi_test(smi, xyzfile, chg, fcmode=True, metal_centers=metal_centers)


def test_6():
    """Metal Center Test 2

    This example shows that the haptic coordination is not supported yet.
    """
    smi = "THIS SHOULD FAIL"
    chg = 0
    metal_centers = [0]

    xyzfile = (
        (pathlib.Path(__file__).parent / "./zirconium_ligand.xyz").absolute().as_posix()
    )  # haptic coordination is not supported.
    assert (
        xyz2smi_test(smi, xyzfile, chg, fcmode=True, metal_centers=metal_centers)
        is False
    )


def test_7():
    """Metal Center Test 3"""
    smi = "[H]c1c(C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])c([H])c2c(c1C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])O[Ti-]13([O+]4C([H])([H])C([H])([H])C([H])([H])C([H])([H])C4([H])[H])Oc4c(C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])c([H])c(C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])c([H])c4C([H])([H])[N+]1(C2([H])[H])C([H])([H])c1c([H])c(C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])c([H])c(C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])[H])c1O3"
    chg = 1  # this turns out to be a cation
    metal_centers = [0]

    xyzfile = (
        (pathlib.Path(__file__).parent / "./titanium_ligand.xyz").absolute().as_posix()
    )
    assert (
        xyz2smi_test(smi, xyzfile, chg, fcmode=True, metal_centers=metal_centers)
        is True
    )


if __name__ == "__main__":

    def run_test(test):
        try:
            test()
        except AssertionError:
            print(f"Test {test.__name__} Failed")
        else:
            print(f"Test {test.__name__} Passed")

    run_test(test_1)
    run_test(test_2)
    run_test(test_3)
    run_test(test_4)
    run_test(test_5)
    run_test(test_6)
    run_test(test_7)
