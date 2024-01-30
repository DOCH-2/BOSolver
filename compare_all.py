from typing import TextIO
from rdkit import Chem
from rdkit.Chem.rdchem import ResonanceMolSupplier
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem import rdmolops
from rdkit.Chem import Draw
from openbabel import pybel, openbabel
import sys

sys.path.append("/home/leejinwon")

from acerxn import chem
from acerxn import process
import numpy as np
import compute_chg_and_bo_gurobi
import time
import subprocess


def contain_atom_not_allowed(rd_mol):
    allowed = ["c", "h", "b", "o", "n", "li", "mg", "si", "p", "s", "f", "na", "co", "rh", "ni", "ti", "fe", "cl", "br", "bb", "lg", "pd", "i"]
    for atom in rd_mol.GetAtoms():
        if not (atom.GetSymbol().lower() in allowed):
            return True
    return False

def restore_rdmol(z_list, bond_dict):
    e_mol = Chem.rdchem.EditableMol()

    # atom construction
    for z in z_list:
        e_mol.AddAtom(Chem.rdchem.Atom(int(z)))
    # bond construction
    for i, order in bond_dict.keys():
        e_mol.AddBond(int(i[0]), int(i[1]), Chem.rdchem.BondType.values[int(order)])

    return e_mol.GetMol()


def get_fail_type(rd_mol, ori_mol):
    bo = Chem.GetAdjacencyMatrix(rd_mol, useBO=True)
    nbonds = int(bo.sum() / 2)
    formal_chgs = np.array([atom.GetFormalCharge() for atom in rd_mol.GetAtoms()])
    chg_sepa = np.sum(np.abs(formal_chgs))

    bo_ori = Chem.GetAdjacencyMatrix(ori_mol, useBO=True)
    nbonds_ori = int(bo_ori.sum() / 2)
    formal_chgs_ori = np.array([atom.GetFormalCharge() for atom in ori_mol.GetAtoms()])
    chg_sepa_ori = np.sum(np.abs(formal_chgs_ori))

    boGE = nbonds >= nbonds_ori
    chgLE = chg_sepa <= chg_sepa_ori

    # Type I: boGE & chgLE
    # Type II: boGE & not chgLE
    # Type III: not boGE & chgLE
    # Type IV: not boGE & not chgLE
    return 4 - (2 ** int(boGE) + int(chgLE))


def restore_rdkit(xyzstr: str, chg: int):
    raw_mol = Chem.MolFromXYZBlock(xyzstr)
    rd_mol = Chem.Mol(raw_mol)
    t = time.time()
    try:
        rdDetermineBonds.DetermineBonds(rd_mol, charge=chg)
        rdkit_time = time.time() - t
        for bond in rd_mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                bond.SetStereo(Chem.BondStereo.STEREONONE)
            elif bond.GetBondType() == Chem.BondType.SINGLE:
                bond.SetBondDir(Chem.BondDir.NONE)
        for atom in rd_mol.GetAtoms():
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        smi_rdkit = Chem.MolToSmiles(rd_mol, canonical=True)
    except:
        rdkit_time = time.time() - t
        smi_rdkit = "BOND ASSIGN FAIL"
        # return smi_rdkit, rdkit_time
    # final check
    try:
        Chem.SanitizeMol(rd_mol)
    except:
        smi_rdkit = "SANITIZE FAIL"
        # return smi_rdkit, rdkit_time

    return smi_rdkit, rdkit_time


def restore_ac2bo(xyzstr: str, chg: int):
    ace_mol = chem.Molecule()
    atom_list = []

    lines = xyzstr.strip().split("\n")
    for i, line in enumerate(lines):
        if i < 2:
            continue
        else:
            content = line.strip()
            atom_line = content.split()
            # atomic_number = int(atom_line[0])
            element_symbol = atom_line[0]
            x = float(atom_line[1])
            y = float(atom_line[2])
            z = float(atom_line[3])
            new_atom = chem.Atom(element_symbol)
            new_atom.x = x
            new_atom.y = y
            new_atom.z = z
            atom_list.append(new_atom)
    ace_mol.atom_list = atom_list
    # At least make adjacency
    ace_mol.adj_matrix = process.get_adj_matrix_from_distance(ace_mol, coeff=1.1)
    ace_mol.atom_feature["chg"] = None
    ace_mol.bo_matrix = None

    t = time.time()
    chg_list, bo_matrix = compute_chg_and_bo_gurobi.compute_chg_and_bo(
        ace_mol, chg, resolve=True
    )
    ac2bo_time = time.time() - t

    if chg_list is None:
        smi_ac2bo = "BOND ASSIGN FAIL"
        return smi_ac2bo, ac2bo_time

    ace_mol.atom_feature["chg"] = chg_list
    ace_mol.bo_matrix = bo_matrix
    rd_actobo = ace_mol.get_rd_mol()
    # rd_actobo = Chem.RemoveHs(rd_actobo, implicitOnly=True)

    sanflg = (
        rdmolops.SanitizeFlags.SANITIZE_CLEANUP
        | rdmolops.SanitizeFlags.SANITIZE_PROPERTIES
        | rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS
        | rdmolops.SanitizeFlags.SANITIZE_KEKULIZE
        | rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS
        | rdmolops.SanitizeFlags.SANITIZE_SETAROMATICITY
        | rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION
        | rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
    )
    # print("sansansan", sanflg)
    try:
        Chem.SanitizeMol(rd_actobo)
        smi_ac2bo = Chem.MolToSmiles(rd_actobo, canonical=True)
    except:
        # smi_ac2bo = Chem.MolToSmiles(rd_actobo)
        # Draw.MolToFile(rd_actobo, "SANITIZE.png")
        smi_ac2bo = "SANITIZE FAIL"
        pass

    return smi_ac2bo, ac2bo_time


def restore_openbabel(xyzstr: str, chg: int):
    ob_mol = pybel.readstring("xyz", xyzstr)
    ob_mol.OBMol.SetTotalCharge(chg)
    ob_mol.OBMol.ConnectTheDots()
    t = time.time()
    try:
        ob_mol.OBMol.PerceiveBondOrders()
        openbabel_time = time.time() - t
    except:
        openbabel_time = time.time() - t
    try:
        smi_tmp = ob_mol.write("smi").strip()
        rd_mol = Chem.MolFromSmiles(smi_tmp)
        Chem.SanitizeMol(rd_mol)
        smi_openbabel = Chem.MolToSmiles(rd_mol)
    except:
        smi_openbabel = smi_tmp
        return smi_openbabel, openbabel_time

    return smi_openbabel, openbabel_time


def restore_indigox(xyzstr: str, chg: int):
    result = subprocess.run(
        ["python", "run_ix.py", xyzstr, str(chg)], capture_output=True, text=True
    )
    if result.returncode != 0:
        indigox_time = 10
        smi_indigox = ""
    else:
        output = result.stdout.strip().split("\n")
        # print(f"indigox:\n{output}")
        count = int(output[0])
        indigox_time = float(output[2])
        smi_indigox = "" if count == 0 else output[4]

    return smi_indigox, indigox_time


def compare_all():
    start, num, set_name = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    num_calculated = 0
    smiles_list = open(f"test_{set_name}.txt", "r")
    ft1 = open(f"{set_name}/ACtoBO_failtype1.txt", "w")
    ft2 = open(f"{set_name}/ACtoBO_failtype2.txt", "w")
    ft3 = open(f"{set_name}/ACtoBO_failtype3.txt", "w")
    ft4 = open(f"{set_name}/ACtoBO_failtype4.txt", "w")
    fts = [ft1, ft2, ft3, ft4]
    ft_bond = open(f"{set_name}/ACtoBO_bondfail.txt", "w")
    ft_sani = open(f"{set_name}/ACtoBO_sanifail.txt", "w")

    rdft1 = open(f"{set_name}/RDKit_failtype1.txt", "w")
    rdft2 = open(f"{set_name}/RDKit_failtype2.txt", "w")
    rdft3 = open(f"{set_name}/RDKit_failtype3.txt", "w")
    rdft4 = open(f"{set_name}/RDKit_failtype4.txt", "w")
    rdfts = [rdft1, rdft2, rdft3, rdft4]
    rdft_bond = open(f"{set_name}/RDKit_bondfail.txt", "w")
    rdft_sani = open(f"{set_name}/RDKit_sanifail.txt", "w")

    RDKit = open(f"{set_name}/RDKit.txt", "w")
    ACtoBO = open(f"{set_name}/ACtoBO.txt", "w")
    # Openbabel = open(f"{set_name}/Openbabel.txt", "w")
    # Indigox = open(f"{set_name}/Indigox.txt", "w")

    line = 0
    ended = 0
    success_rdkit = 0
    success_actobo = 0
    # success_openbabel = 0
    # success_indigox = 0

    total_rdkit = 0
    total_actobo = 0
    # total_openbabel = 0
    # total_indigox = 0

    rdkit_list = []
    actobo_list = []
    # openbabel_list = []
    # indigox_list = []

    for smiles in smiles_list:
        line += 1
        if line < start or len(smiles) > 160:
            continue
        if ended == num:
            break
        print("Current line:", line)
        print("smi:", smiles.strip())
        print()

        try:
            molecule = chem.Molecule(smiles)
        except:
            continue
        conformers = molecule.sample_conformers(1)
        if len(conformers) == 0:
            continue
        rd_org = molecule.get_rd_mol()
        
        # use ONLY organic compunds
        if contain_atom_not_allowed(rd_org):
            continue

        # use ONLY Sanitizable SMILES
        try:
            Chem.SanitizeMol(rd_org)
        except:
            continue
        smi_original = Chem.MolToSmiles(rd_org, canonical=True)
        num_calculated += 1

        try:
            molecule = chem.Molecule(smi_original)
        except:
            continue
        conformers = molecule.sample_conformers(1)
        if len(conformers) == 0:
            continue
        conformer = conformers[0]

        xyzstr = conformer.get_content()
        chg = molecule.get_chg()
        # conformer.write_geometry("input.xyz")

        ### rdkit ###
        smi_rdkit, rdkit_time = restore_rdkit(xyzstr, chg)

        ### ACtoBO ###
        smi_actobo, actobo_time = restore_ac2bo(xyzstr, chg)

        ### openbabel ###
        # smi_openbabel, openbabel_time = restore_openbabel(xyzstr, chg)

        ### indigox ###
        # smi_indigox, indigox_time = restore_indigox(xyzstr, chg)

        ### Comparison ###
        flag_rdkit = smi_original == smi_rdkit
        flag_actobo = smi_original == smi_actobo
        # flag_openbabel = smi_original == smi_openbabel
        # flag_indigox = smi_original == smi_indigox
        res_flag = (
            Chem.KEKULE_ALL
            | Chem.ALLOW_INCOMPLETE_OCTETS
            | Chem.UNCONSTRAINED_CATIONS
            | Chem.UNCONSTRAINED_CATIONS
            | Chem.ALLOW_CHARGE_SEPARATION
        )
        suppl = ResonanceMolSupplier(rd_org, flags=res_flag)
        for resMol in suppl:
            smi_res = Chem.MolToSmiles(resMol)
            if smi_res == smi_rdkit:
                flag_rdkit = True
            if smi_res == smi_actobo:
                flag_actobo = True
            # if smi_res == smi_openbabel:
            #    flag_openbabel = True
            # if smi_res == smi_indigox:
            #    flag_indigox = True

        if flag_rdkit:
            success_rdkit += 1
        else:
            rdkit_list.append((smi_original, smi_rdkit))
            RDKit.write("Original: " + smi_original + "\n")
            RDKit.write("RDKit: " + smi_rdkit + "\n\n")

            if smi_rdkit == "BOND ASSIGN FAIL":
                rdft_bond.write("Original: " + smi_original + "\n")
                rdft_bond.write("RDKit: " + smi_rdkit + "\n\n")
            elif smi_rdkit == "SANITIZE FAIL":
                rdft_sani.write("Original: " + smi_original + "\n")
                rdft_sani.write("RDKit: " + smi_rdkit + "\n\n")
            else:
                failtype = get_fail_type(
                    Chem.MolFromSmiles(smi_rdkit), Chem.MolFromSmiles(smi_original)
                )
                rdfts[failtype - 1].write("Original: " + smi_original + "\n")
                rdfts[failtype - 1].write("RDKit: " + smi_rdkit + "\n\n")

        if flag_actobo:
            success_actobo += 1
        else:
            actobo_list.append((smi_original, smi_actobo))
            ACtoBO.write("Original: " + smi_original + "\n")
            ACtoBO.write("ACtoBO: " + smi_actobo + "\n\n")

            if smi_actobo == "BOND ASSIGN FAIL":
                ft_bond.write("Original: " + smi_original + "\n")
                ft_bond.write("ACtoBO: " + smi_actobo + "\n\n")
            elif smi_actobo == "SANITIZE FAIL":
                ft_sani.write("Original: " + smi_original + "\n")
                ft_sani.write("ACtoBO: " + smi_actobo + "\n\n")
            else:
                failtype = get_fail_type(
                    Chem.MolFromSmiles(smi_actobo), Chem.MolFromSmiles(smi_original)
                )
                fts[failtype - 1].write("Original: " + smi_original + "\n")
                fts[failtype - 1].write("ACtoBO: " + smi_actobo + "\n\n")

        # if flag_openbabel:
        #    success_openbabel += 1
        # else:
        #    openbabel_list.append((smi_original, smi_openbabel))
        #    Openbabel.write("Original: " + smi_original + "\n")
        #    Openbabel.write("Openbabel: " + smi_openbabel + "\n\n")

        # if flag_indigox:
        # success_indigox += 1
        # else:
        # indigox_list.append((smi_original, smi_indigox))
        # Indigox.write("Original: " + smi_original + "\n")
        # Indigox.write("Indigox: " + smi_indigox + "\n\n")

        ended += 1
        print(str(ended) + "/" + str(num))

        color = "\033[92m" if flag_rdkit else "\033[91m"
        print(color + "RDKit time: " + str(rdkit_time))

        color = "\033[92m" if flag_actobo else "\033[91m"
        print(color + "ACtoBO time: " + str(actobo_time))

        # color = "\033[92m" if flag_openbabel else "\033[91m"
        # print(color + "Openbabel time: " + str(openbabel_time))

        # color = "\033[92m" if flag_indigox else "\033[91m"
        # print(color + "indigox time: " + str(indigox_time))

        print("\033[0m")

        total_rdkit += rdkit_time
        total_actobo += actobo_time
        # total_openbabel += openbabel_time
        # total_indigox += indigox_time

    print(
        "RDKit"
        + " (success/total):"
        + f" {success_rdkit}/{num}"
        + f" ({100 * success_rdkit / num:.1f}%)"
        + " (success/calculated): "
        + f" {success_rdkit}/{num_calculated}"
        + f" ({100 * success_rdkit / num_calculated:.1f}%)"
        + f" average TIME"
        + f"{total_rdkit / num}"
    )
    print(
        "AC2BO"
        + " (success/total):"
        + f" {success_actobo}/{num}"
        + f" ({100 * success_actobo/ num:.1f}%)"
        + " (success/calculated): "
        + f" {success_actobo}/{num_calculated}"
        + f" ({100 * success_actobo/ num_calculated:.1f}%)"
        + f" average TIME"
        + f" {total_actobo / num}"
    )
    # print(
    # "Openbabel: "
    # + str(success_openbabel)
    # + "/"
    # + str(num)
    # + " ("
    # + str(100 * success_openbabel / num)
    # + "%) "
    # + str(total_openbabel / num)
    # )
    # print(
    # "indigox: "
    # + str(success_indigox)
    # + "/"
    # + str(num)
    # + " ("
    # + str(100 * success_indigox / num)
    # + "%) "
    # + str(total_indigox / num)
    # )
    print()

    smiles_list.close()
    # bo_notsmall_chg_notbig_list.close()
    # bo_notsmall_chg_big_list.close()
    # bo_small_chg_notbig_list.close()
    # bo_small_chg_big_list.close()
    # rdkit_bo_notsmall_chg_notbig_list.close()
    # rdkit_bo_notsmall_chg_big_list.close()
    # rdkit_bo_small_chg_notbig_list.close()
    # rdkit_bo_small_chg_big_list.close()

    RDKit.close()
    ACtoBO.close()
    # Openbabel.close()
    # Indigox.close

    # flag = input("Print Counterexamples (Y/N) ")
    flag = "Y"
    print("Print Counterexamples")
    print()
    if flag == "Y":
        print("===RDKit===")
        for s1, s2 in rdkit_list:
            print("Original: ", s1)
            print("RDKit: ", s2)
            print()
        print("===ACtoBO===")
        for s1, s2 in actobo_list:
            print("Original: ", s1)
            print("ACtoBO: ", s2)
            print()
        # print("===Openbabel===")
        # for s1, s2 in openbabel_list:
        # print("Original: ", s1)
        # print("Openbabel: ", s2)
        # print()
        # print("===Indigox===")
        # for s1, s2 in indigox_list:
        # print("Original: ", s1)
        # print("Indigox: ", s2)
        # print()


if __name__ == "__main__":
    compare_all()
