"""parse raw data (json file) into individual xyz files"""

import ijson
import numpy as np
import pandas as pd

from rdkit import Chem
from openbabel import openbabel as ob

PT = Chem.GetPeriodicTable()

with open("./1-298189.json", "r") as f:
    prefix1 = "cid"
    prefix2 = "coordinates.item"
    prefix3 = "atomic-numbers.item"
    prefix4 = "number-of-atoms"

    prefix5 = "pubchem-inchi"
    prefix6 = "pubchem-isomeric-smiles"
    prefix7 = "pubchem-obabel-canonical-smiles"
    prefix8 = "pubchem-charge"

    prefix9 = "pm6-obabel-canonical-smiles"
    prefix10 = "obabel-inchi"
    prefix11 = "charge"

    parser = ijson.parse(f)

    cid = 0
    coordinates = []
    z_list = []
    nAtoms = 0
    chg = 0
    smi_ori = ""
    smi_obb = ""
    pubchem_inchi = ""
    smi_pubchem = ""

    i = 0
    for prefix, event, value in parser:
        if prefix.endswith(prefix1):
            # check if it's really neutral...
            if cid > 0:
                conv = ob.OBConversion()
                conv.SetInFormat("smi")
                obmol = ob.OBMol()
                conv.ReadString(obmol, smi_ori)
                chg_from_smi = obmol.GetTotalCharge()

                if chg_from_smi != chg:
                    print(
                        f"CID {cid}: Wrong chg were used. overwrites charge to charge from SMILES"
                    )
                    chg = chg_from_smi
                    break

            # flush cached data
            if chg == 0:
                category = "neutral"
            elif chg > 0:
                category = "cation"
            else:
                category = "anion"
            coordinates = np.array(coordinates).reshape(-1, 3)
            with open(f"{category}/{cid}.xyz", "w") as f:
                f.write(str(nAtoms) + "\n")
                f.write(f"chg:{chg} smi_ori:{smi_ori} smi_obb:{smi_obb}\n")
                for z, coord in zip(z_list, coordinates):
                    f.write(
                        f"{PT.GetElementSymbol(int(z))} {' '.join([str(x) for x in coord.tolist()])}"
                    )
                    f.write("\n")

            # reset cache
            cid = 0
            coordinates = []
            z_list = []
            nAtoms = 0
            chg = 0
            smi_ori = ""
            smi_obb = ""
            pubchem_inchi = ""
            smi_pubchem = ""

            # new data
            cid = value

        elif prefix.endswith(prefix2):
            coordinates.append(value)
        elif prefix.endswith(prefix3):
            z_list.append(value)
        elif prefix.endswith(prefix4):
            nAtoms = value
        elif prefix.endswith(prefix5):
            pubchem_inchi = value
        elif prefix.endswith(prefix6):
            smi_pubchem = value
        elif prefix.endswith(prefix8):
            chg = value
        elif prefix.endswith(prefix7):
            smi_ori = value
        elif prefix.endswith(prefix9):
            smi_obb = value
