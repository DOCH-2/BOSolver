import os.path as osp
import argparse
import subprocess

from convert import (
    obabel_charge,
    indigox_charge,
    xyz2mol_charge,
    main_charge,
    obabel_nocharge,
    xyz2mol_nocharge,
    indigox_nocharge,
    main_nocharge,
    nocharge_chgs,
)
from grade_noAns import check_All


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

name = osp.basename(xyzfile)

if chg is not None:
    if useOBABEL:
        f = obabel_charge
    elif useXYZ2MOL:
        f = xyz2mol_charge
    elif useINDIGOX:
        f = indigox_charge
    else:
        f = main_charge
    mol, t = f(xyzfile, chg)
    # name, chg, chgConserve, noRadical, sanitize, boSum, chgSepa, nAromatic, nOctet, time
    print(name, chg, *check_All(mol, chg), f"{t:.5f}", sep="\t")

else:
    if useINDIGOX:
        """indigox cannot be called sequentially.
        It will raise an segmentation fault if the process is not closed after bond order assignment.
        It seems that memory leak is involved in the process.
        Therefore, instead of calling indigox multiple times in the same process,
        the script will call indigox for each charge state of the molecule in different subprocesses.
        """
        chgs = nocharge_chgs(xyzfile)
        args = [
            [
                "python",
                __file__,
                xyzfile,
                str(chg),
                "--indigo",
            ]
            for chg in chgs
        ]
        outputs = []
        for arg in args:
            out = subprocess.run(arg, capture_output=True, text=True).stdout
            outputs.append(out)

        for out in outputs:
            print(out)

    else:
        if useOBABEL:
            f = obabel_nocharge
        elif useXYZ2MOL:
            f = xyz2mol_nocharge
        else:
            f = main_nocharge
        results = f(xyzfile)
        for result in results:
            (mol, t), chg = result
            # name, chg, noRadical, sanitize, boSum, chgSepa, nAromatic, time
            print(name, chg, *check_All(mol, chg), f"{t:.5f}", sep="\t")
