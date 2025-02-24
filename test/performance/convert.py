import os
import os.path as osp
import sys
import argparse
import time

import numpy as np

from rdkit import Chem

# Algorithms
from openbabel import openbabel as ob
import indigox as ix
import xyz2mol

from BOSolver.bosolve import perceiveConn, assignBO

from BOSolver.utils.coord2adj import CovalentRadius

ADJ_COEFF = 1.15


def obabel_charge(xyzfile, chg=None):
    # chg is never used in this function
    name = osp.basename(xyzfile)
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("xyz", "smi")
    conv.SetOptions("c", conv.OUTOPTIONS)
    conv.SetOptions("h", conv.OUTOPTIONS)
    conv.SetOptions("i", conv.OUTOPTIONS)

    obmol = ob.OBMol()
    conv.ReadFile(obmol, xyzfile)

    obmol.ConnectTheDots()

    # Algorithm: Openbabel #########################
    try:
        s = time.time()
        obmol.PerceiveBondOrders()
        e = time.time()
    except Exception as exc:
        t = -1.0
        rdmol = None
        print(f"Error ({name}, {chg}): {exc} ({type(exc)=})", file=sys.stderr)
    else:
    ##############################################
        t = e - s
        _rdmol = Chem.RWMol()
        for obatom in ob.OBMolAtomIter(obmol):
            _rdatom = Chem.Atom(obatom.GetAtomicNum())
            _rdatom.SetFormalCharge(obatom.GetFormalCharge())
            _rdmol.AddAtom(_rdatom)
        for obbond in ob.OBMolBondIter(obmol):
            bo = obbond.GetBondOrder()
            if bo == 1:
                bo = Chem.BondType.SINGLE
            elif bo == 2:
                bo = Chem.BondType.DOUBLE
            elif bo == 3:
                bo = Chem.BondType.TRIPLE
            elif bo == 1.5:
                bo = Chem.BondType.AROMATIC
            else:
                bo = Chem.BondType.UNSPECIFIED
            _rdmol.AddBond(
                obbond.GetBeginAtomIdx() - 1,
                obbond.GetEndAtomIdx() - 1,
                bo,
            )
        rdmol = _rdmol.GetMol()
    return rdmol, t


def xyz2mol_charge(xyzfile, chg):
    name = osp.basename(xyzfile)
    atoms, _, xyz_coords = xyz2mol.read_xyz_file(xyzfile)

    # Algorithm: xyz2mol #########################
    try:
        s = time.time()
        # sometimes xyz2mol prints error code to stdout
        with open(os.devnull, 'w') as null:
            sys.stdout = null
            mols = xyz2mol.xyz2mol(atoms, xyz_coords, chg)
        sys.stdout = sys.__stdout__
        e = time.time()
        t = e - s
        mol = mols[0]
    except IndexError:
        mol = None
        t = -1.0
        print(f"Error ({name}, {chg}): conversion failed", file=sys.stderr)
    except Exception as exc:
        mol = None
        t = -1.0
        print(f"Error ({name}, {chg}): {exc} ({type(exc)=})", file=sys.stderr)
    ##############################################
    return mol, t


def indigox_charge(xyzfile, chg):
    pt = ix.PeriodicTable()

    name = osp.basename(xyzfile)
    mol = Chem.MolFromXYZFile(xyzfile)
    perceive_algo = CovalentRadius(relTol=ADJ_COEFF, absTol=0.0)
    mol = perceive_algo(mol, algorithm=perceive_algo)

    ixmol = ix.Molecule()

    for i, atom in enumerate(mol.GetAtoms()):
        atom = ixmol.NewAtom(i, pt.GetElement(atom.GetAtomicNum()))
        atom.SetName(f"{i}")

    adj = Chem.GetAdjacencyMatrix(mol)
    for i, j in np.stack(np.where(adj > 0)).T:
        if i < j:
            ixmol.NewBond(ixmol.GetAtomIndex(i), ixmol.GetAtomIndex(j))

    opts = ix.Options.AssignElectrons
    opts.ALGORITHM = opts.Algorithm.FPT
    opts.FPT.ADD_EDGES_TO_TD = False
    opts.FPT.MINIMUM_PROPAGATION_DEPTH = 1
    opts.USE_ELECTRON_PAIRS = False

    ixmol.SetTotalCharge(chg)

    # Algorithm: IndigoX #########################
    try:
        s = time.time()
        counts = ixmol.AssignElectrons()
        e = time.time()
        assert counts > 0
    except AssertionError:
        t = -1.0
        rdmol = None
        print(f"Error ({name}, {chg}): conversion failed", file=sys.stderr)
    except Exception as exc:
        t = -1.0
        rdmol = None
        print(f"Error ({name}, {chg}): {exc} ({type(exc)=})", file=sys.stderr)
    ##############################################
    else:
        t = e - s
        _rdmol = Chem.RWMol()
        ixmol.ApplyElectronAssignment(0)

        for ixatom in ixmol.GetAtoms():
            _rdatom = Chem.Atom(ixatom.GetElement().GetAtomicNumber())
            _rdatom.SetFormalCharge(ixatom.GetFormalCharge())
            _rdmol.AddAtom(_rdatom)
        for ixbond in ixmol.GetBonds():
            ixbo = ixbond.GetOrder()
            if ixbo == ix.BondOrder.SINGLE_BOND:
                bo = Chem.BondType.SINGLE
            elif ixbo == ix.BondOrder.DOUBLE_BOND:
                bo = Chem.BondType.DOUBLE
            elif ixbo == ix.BondOrder.TRIPLE_BOND:
                bo = Chem.BondType.TRIPLE
            elif ixbo == ix.BondOrder.AROMATIC_BOND:
                bo = Chem.BondType.AROMATIC
            else:
                bo = Chem.BondType.UNSPECIFIED

            _rdmol.AddBond(
                ixbond.GetSourceAtom().GetIndex(),
                ixbond.GetTargetAtom().GetIndex(),
                bo,
            )

        rdmol = _rdmol.GetMol()


    return rdmol, t


def main_charge(xyzfile, chg):
    name = osp.basename(xyzfile)
    mol = Chem.MolFromXYZFile(xyzfile)
    perceive_algo = CovalentRadius(relTol=ADJ_COEFF, absTol=0.0)
    mol = perceiveConn(mol, algorithm=perceive_algo)

    # Algorithm: BOSolver #########################
    try:
        s = time.time()
        mol = assignBO(mol, chg)
        e = time.time()
    except AssertionError:
        t = -1.0
        # smi = None
        mol = None
        print(f"Error ({name}, {chg}): conversion failed", file=sys.stderr)
    except Exception as exc:
        t = -1.0
        # smi = None
        mol = None
        print(f"Error ({name}, {chg}): {exc} ({type(exc)=})", file=sys.stderr)
    else:
    ###############################################
        t = e - s
    return mol, t


def obabel_nocharge(xyzfile):
    mol, t = obabel_charge(xyzfile, chg=None)
    chg = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    return [((mol, t), chg)]


def xyz2mol_nocharge(xyzfile):
    atoms, _, xyz_coords = xyz2mol.read_xyz_file(xyzfile)

    z_sum = sum(atoms)
    if z_sum % 2 == 0:
        chgs = [0, 2, -2]
    else:
        chgs = [1, -1]
    return [(xyz2mol_charge(xyzfile, chg), chg) for chg in chgs]


def indigox_nocharge(xyzfile):
    atoms, _, xyz_coords = xyz2mol.read_xyz_file(xyzfile)

    z_sum = sum(atoms)
    if z_sum % 2 == 0:
        chgs = [0, 2, -2]
    else:
        chgs = [1, -1]
    return [(indigox_charge(xyzfile, chg), chg) for chg in chgs]


def main_nocharge(xyzfile):
    atoms, _, xyz_coords = xyz2mol.read_xyz_file(xyzfile)

    z_sum = sum(atoms)
    if z_sum % 2 == 0:
        chgs = [0, 2, -2]
    else:
        chgs = [1, -1]
    return [(main_charge(xyzfile, chg), chg) for chg in chgs]


def nocharge_chgs(xyzfile):
    atoms, _, xyz_coords = xyz2mol.read_xyz_file(xyzfile)

    z_sum = sum(atoms)
    if z_sum % 2 == 0:
        chgs = [0, 2, -2]
    else:
        chgs = [1, -1]
    return chgs


if __name__ == "__main__":
    import os.path as osp
    from rdkit.Chem import Draw

    p = argparse.ArgumentParser()
    p.add_argument("xyzfile", type=str)
    p.add_argument("chg", type=int, default=0)
    conv_method = p.add_mutually_exclusive_group()
    conv_method.add_argument("--bosolver", dest="bosolver", action="store_true")
    conv_method.add_argument("--obabel", dest="obabel", action="store_true")
    conv_method.add_argument("--indigo", dest="indigo", action="store_true")
    conv_method.add_argument("--xyz2mol", dest="xyz2mol", action="store_true")

    args = p.parse_args()
    xyzfile = args.xyzfile
    chg = args.chg
    useBOSOLV = args.bosolver
    useOBABEL = args.obabel
    useINDIGOX = args.indigo
    useXYZ2MOL = args.xyz2mol

    if useOBABEL:
        mol, _ = obabel_charge(xyzfile, chg)
    elif useXYZ2MOL:
        mol, _ = xyz2mol_charge(xyzfile, chg)
        Chem.AllChem.Compute2DCoords(mol)
    elif useINDIGOX:
        mol, _ = indigox_charge(xyzfile, chg)
    else:
        mol, _ = main_charge(xyzfile, chg)

    Chem.SanitizeMol(mol)
    for atom in mol.GetAtoms():
        if atom.GetTotalNumHs() > 0:
            atom.SetNumRadicalElectrons(atom.GetTotalNumHs())
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(True)
    print(Chem.MolToSmiles(mol, isomericSmiles=False))