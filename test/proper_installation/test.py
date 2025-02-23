from rdkit import Chem
from rdkit.Chem import AllChem
import sys

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

    xyzfile = "./platinum_ligand.xyz"
    assert xyz2smi_test(smi, xyzfile, chg, fcmode=True, metal_centers=metal_centers)


def test_6():
    """Metal Center Test 2"""
    smi = "THIS SHOULD FAIL"
    chg = 0
    metal_centers = [0]

    xyzfile = "./zirconium_ligand.xyz"  # haptic coordination is not supported.
    assert (
        xyz2smi_test(smi, xyzfile, chg, fcmode=True, metal_centers=metal_centers)
        is False
    )


# FIXME: test for titanium_ligand.xyz is not working
# BOSolver cannot solve a system with odd number of electrons
def test_7():
    """Odd number electrons Test"""
    smi = "THIS FAILS"
    chg = 0
    metal_centers = [0]

    xyzfile = "./titanium_ligand.xyz"
    try:
        xyz2smi_test(smi, xyzfile, chg, fcmode=True, metal_centers=metal_centers)
    except RuntimeError:
        # it fails...
        assert True
    else:
        assert False


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
