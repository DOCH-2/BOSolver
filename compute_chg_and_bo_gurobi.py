from typing import List, Dict, Tuple
from itertools import combinations as comb

from rdkit import Chem

from gurobipy import GRB, Model, LinExpr, Env
from acerxn import chem
import numpy as np


def get_ring_info(z_list, adj_matrix):
    new_z = list(np.where(z_list > 1, 6, 1))

    # print(new_adj)
    chg_list = np.zeros(len(new_z))
    new = chem.Molecule([new_z, adj_matrix, None, chg_list])
    new_rd = new.get_rd_mol()
    Chem.SanitizeMol(new_rd)

    sssrs = Chem.GetSymmSSSR(new_rd)
    atoms_in_ring = [list(sssrs[i]) for i in range(len(sssrs))]
    bond_rings = new_rd.GetRingInfo().BondRings()
    bonds_in_ring = [
        list(
            map(
                lambda x: (
                    new_rd.GetBondWithIdx(x).GetBeginAtomIdx(),
                    new_rd.GetBondWithIdx(x).GetEndAtomIdx(),
                ),
                ring,
            )
        )
        for ring in bond_rings
    ]
    # print("bondsInRing:", bondsInRing)
    # print(len(bondsInRing))
    # print("AtomInRings:", [list(sssrs[i]) for i in range(len(sssrs))])
    # print(len([list(sssrs[i]) for i in range(len(sssrs))]))

    return atoms_in_ring, bonds_in_ring


def get_lists(molecule: chem.Molecule):
    # period, group, adj
    period_list, group_list = molecule.get_period_group_list()
    adj_matrix = np.copy(molecule.get_matrix("adj"))
    adj_list = np.sum(adj_matrix, axis=1)

    # neighbor, bond, bond mapping
    neighbor_list = molecule.get_neighbor_list()
    bond_list = molecule.get_bond_list(False)
    bond_mapping = {key: val for val, key in enumerate(bond_list)}

    # valence, atomic number
    ve_list = np.zeros_like(group_list)
    z_list = molecule.get_z_list()
    for i in range(len(group_list)):
        if period_list[i] == 1:
            ve_list[i] = 2
        else:
            # not considering expanded octet here
            ve_list[i] = 8

    # ring membership
    ring_list, ring_bond_list = get_ring_info(z_list, adj_matrix)

    return (
        period_list,
        group_list,
        z_list,
        ve_list,
        adj_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        ring_list,
        ring_bond_list,
    )


def get_expanded_list(period_list, ve_list, chg_list, ring_list):
    ring_members = np.unique(sum(ring_list, []))

    in_ring = np.zeros_like(period_list)
    if len(ring_members) > 0:
        in_ring[ring_members] = 1
    in_ring = in_ring.astype(bool)

    expanded_idx = (period_list > 2) & ~in_ring
    eve_list = np.copy(ve_list)
    eve_list[expanded_idx] += 2 * np.where(chg_list > 0, chg_list, 0)[expanded_idx]

    return eve_list

def get_modified_list(period_list, eve_list, chg_list, bo_dict, ring_list, ring_bond_list):
    mve_list = np.copy(eve_list) 
    ring_members = np.unique(sum(ring_list, []))
    in_ring = np.zeros_like(period_list)
    if len(ring_members) > 0:
        in_ring[ring_members] = 1
    in_ring = in_ring.astype(bool)
    
    # CleanUp valence expansion
    # subject:
    # period>2, ring member, one or no double/triple bonds with
    # other ring members
    rbtbc = np.zeros_like(period_list)
    for ring_bonds in ring_bond_list:
        for bond in ring_bonds:
            if bo_dict[bond] > 1:
                rbtbc[bond[0]] += 1
                rbtbc[bond[1]] += 1
    cleanUp_idx = (np.array(period_list) > 2) & in_ring & (rbtbc < 2)
    
    mve_list[cleanUp_idx] += 2 * np.where(chg_list > 0, chg_list, 0)[cleanUp_idx]
    
    return mve_list
    


def maximize_bo(
    atom_num,
    bond_num,
    group_list,
    bond_list,
    bond_mapping,
    ve_list,
    neighbor_list,
    chg_mol,
    eIsEven,
    **kwargs
):
    # early stop
    if atom_num == 1:
        return np.array([chg_mol]), {}

    ### model construction
    env = Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.start()
    model = Model("maximize_bo", env=env)

    # bo: bond order
    bo = model.addVars(bond_num, lb=1, name="bo", vtype=GRB.INTEGER)

    # t1: formal charge
    # t1[2i]: fc+ | t1[2i+1]: fc-
    # t1[2i] - t1[2i+1] : fc of atom i
    # t1[2i] + t1[2i+1] : abs(fc) of atom i
    t1 = model.addVars(2 * atom_num, lb=0, name="t1", vtype=GRB.INTEGER)

    # t2: formal charge for weighted objective function
    # weight considering electronegativity
    # t2 = model.addVars(2 * atom_num, name="t2", vtype=GRB.CONTINUOUS)

    # even: dummy variable to force no. of electrons even
    even = model.addVars(atom_num, name="even", vtype=GRB.INTEGER)

    ### constraints construction
    chg_constr = LinExpr()  # charge conservation rule
    for i in range(atom_num):
        lp_constr = LinExpr()  # lone pair rule
        ve_constr = LinExpr()  # valence rule

        chg_constr.add(t1[2 * i] - t1[2 * i + 1])
        lp_constr.add(t1[2 * i] - t1[2 * i + 1])
        ve_constr.add(-t1[2 * i] + t1[2 * i + 1])

        # summation over bond
        for j in neighbor_list[i]:
            a, b = i, j
            if a > b:
                a, b = b, a
            lp_constr.add(bo[bond_mapping[(a, b)]])
            ve_constr.add(bo[bond_mapping[(a, b)]])

        model.addConstr(lp_constr <= group_list[i], name=f"lp_{i}")

        # the number of valence electron is less than maximum valence
        model.addConstr(ve_constr + group_list[i] <= ve_list[i], name=f"ve_{i}")

        if eIsEven:
            # the number of valence electron is even number (no radical rule!)
            model.addConstr(ve_constr + group_list[i] == 2 * even[i], name=f"noRad_{i}")
    model.addConstr(chg_constr == chg_mol, name="chg_consv")

    ### optimization
    max_bo_obj = LinExpr()  # bond maximization
    min_fc_obj = LinExpr()  # formal charge minimization
    # min_wfc_obj = LinExpr()  # formal charge minimization (weighted)

    # bond maximization
    for i in range(bond_num):
        max_bo_obj.add(bo[i])
    # (weighted) formal charge minimization
    for i in range(atom_num):
        min_fc_obj.add(t1[2 * i] + t1[2 * i + 1])
        # min_wfc_obj.add((0.1 * group_list[i] + 0.6) * (t2[2 * i] + t2[2 * i + 1]))
        
    bo_priority = 2
    chg_priority = 1
    
    if kwargs.get("mode", "") == "fc":
        bo_priority, chg_priority = chg_priority, bo_priority

    model.setObjectiveN(max_bo_obj, 1, priority=bo_priority, weight=GRB.MAXIMIZE, name="max_bo")
    model.setObjectiveN(min_fc_obj, 2, priority=chg_priority, weight=GRB.MINIMIZE, name="min_fc")
    
    #for test only
    #model.setObjectiveN(max_bo_obj, 0, priority=1, weight=GRB.MAXIMIZE, name="max_bo")
    #model.setObjectiveN(min_fc_obj, 1, priority=2, weight=GRB.MINIMIZE, name="min_fc")
    # if mode == "heuristics":
    #    model.setObjectiveN(min_wfc_obj, 2, 0, GRB.MINIMIZE)

    # Gurobi optimization
    model.setParam(GRB.Param.OutputFlag, 0)
    model.setParam(GRB.Param.TimeLimit, 1)
    #model.write("record.lp")

    model.optimize()

    # error handling
    if model.status != GRB.Status.OPTIMAL:
        return np.zeros(atom_num), {}

    # result record
    nSol = model.SolCount
    # print("nSol", nSol)
    for s in range(nSol):
        model.params.SolutionNumber = s
        #model.write(f"output{s}.sol")
    
    # retrieval
    bo_dict = {}
    chg_list = np.zeros(atom_num, dtype=np.int64)
    for i in range(bond_num):
        bo_dict[bond_list[i]] = int(bo[i].X)
    for i in range(atom_num):
        chg_list[i] = int(t1[2 * i].X) - int(t1[2 * i + 1].X)

    model.close()
    env.close()
    
    return chg_list, bo_dict


def resolve_chg(
    atom_num,
    bond_num,
    group_list,
    bond_list,
    bond_mapping,
    eve_list,
    neighbor_list,
    chg_mol,
    eIsEven,
    alreadyOctet,
    stepIdx=0,
):
    if atom_num == 1:
        return np.array([chg_mol]), {}
    
    ### model construction
    env = Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.start()

    model = Model(f"resolve_chg{stepIdx}", env=env)

    # bo: bond order
    bo = model.addVars(bond_num, lb=1, name="bo", vtype=GRB.INTEGER)

    # t1: formal charge
    t1 = model.addVars(2 * atom_num, lb=0, name="t1", vtype=GRB.INTEGER)
    # t2: formal charge for weighted objective function
    # weight considering electronegativity
    # t2 = model.addVars(2 * atom_num, name="t2", vtype=GRB.CONTINUOUS)

    # even: dummy variable to force no. of electrons even
    even = model.addVars(atom_num, name="even", vtype=GRB.INTEGER)

    ### constraints construction
    chg_constr = LinExpr()  # charge conservation rule
    octet_constr = LinExpr()  #

    for i in range(atom_num):
        lp_constr = LinExpr()  # lone pair rule
        eve_constr = LinExpr()  # expanded valence rule

        chg_constr.add(t1[2 * i] - t1[2 * i + 1])
        lp_constr.add(t1[2 * i] - t1[2 * i + 1])
        eve_constr.add(-t1[2 * i] + t1[2 * i + 1])

        # summation over bond
        for j in neighbor_list[i]:
            a, b = i, j
            if a > b:
                a, b = b, a
            lp_constr.add(bo[bond_mapping[(a, b)]])
            eve_constr.add(bo[bond_mapping[(a, b)]])

        model.addConstr(lp_constr <= group_list[i], name=f"lp_{i}")
        if bool(alreadyOctet[i]):
            model.addConstr(eve_constr + group_list[i] == 8, name=f"eve_{i}")
        else:
            model.addConstr(eve_constr <= eve_list[i] - group_list[i], name=f"eve_{i}")

        if eIsEven:
            # the number of valence electron is even number (no radical rule!)
            model.addConstr(
                eve_constr + group_list[i] == 2 * even[i], name=f"noRad_{i}"
            )

    model.addConstr(chg_constr == chg_mol, name="chg_consv")

    ### optimization
    obj = LinExpr()  # charge separation min & bond max

    # bond maximization
    for i in range(bond_num):
        obj.add(-2 * bo[i])
    # formal charge separation minimization
    for i in range(atom_num):
        obj.add(3 * t1[2 * i] + 3 * t1[2 * i + 1])

    model.setObjective(obj, GRB.MINIMIZE)

    # Gurobi optimization
    model.setParam(GRB.Param.OutputFlag, 0)
    model.setParam(GRB.Param.TimeLimit, 1)
    #model.write(f"record_chg{stepIdx}.lp")

    model.optimize()

    # error handling
    if model.status != GRB.Status.OPTIMAL:
        return np.zeros(atom_num), {}

    # result record
    nSol = model.SolCount
    # print("nSol", nSol)
    for s in range(nSol):
        model.params.SolutionNumber = s
        #model.write(f"output_chg{stepIdx}_{s}.sol")
    
    # retrieval
    bo_dict = {}
    chg_list = np.zeros(atom_num, dtype=np.int64)
    for i in range(bond_num):
        bo_dict[bond_list[i]] = int(bo[i].X)
    for i in range(atom_num):
        chg_list[i] = int(t1[2 * i].X) - int(t1[2 * i + 1].X)

    model.close()
    env.close()

    return chg_list, bo_dict


def compute_chg_and_bo_debug(molecule, chg_mol, resolve=True, cleanUp=True, **kwargs):
    (
        period_list,
        group_list,
        z_list,
        ve_list,
        adj_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        ring_list,
        ring_bond_list,
    ) = get_lists(molecule)
    atom_num, bond_num = len(z_list), len(bond_list)
    eIsEven = int(np.sum(z_list) - chg_mol) % 2 == 0
    resolve_step = 0

    chg_list, bo_dict = maximize_bo(
        atom_num,
        bond_num,
        group_list,
        bond_list,
        bond_mapping,
        ve_list,
        neighbor_list,
        chg_mol,
        eIsEven,
        **kwargs
    )
    # early stop
    if len(bo_dict) == 0:
        return None, None, None, None

    bo_matrix = np.zeros((atom_num, atom_num))
    for p, q in bo_dict.keys():
        bo_matrix[p][q] = bo_dict[(p, q)]
        bo_matrix[q][p] = bo_dict[(p, q)]

    # check charge separation
    chg_sep = np.any(chg_list > 0) and np.any(chg_list < 0)

    # charge resolution
    if resolve and chg_sep:
        bo_sum = np.zeros(atom_num)
        for p, q in bo_dict.keys():
            bo_sum[p] += bo_dict[(p, q)]
            bo_sum[q] += bo_dict[(p, q)]
        alreadyOctet = (group_list + bo_sum - chg_list == 8) & (period_list == 2)
        print("alreadyOctet", np.nonzero(alreadyOctet))

        # ve_list is now expanded Valence list
        eve_list = get_expanded_list(period_list, ve_list, chg_list, ring_list)
        new_chg_list, new_bo_dict = resolve_chg(
            atom_num,
            bond_num,
            group_list,
            bond_list,
            bond_mapping,
            eve_list,
            neighbor_list,
            chg_mol,
            eIsEven,
            alreadyOctet,
            resolve_step,
        )
        resolve_step += 1
        
        if cleanUp:
            mve_list = get_modified_list(period_list, eve_list, new_chg_list, new_bo_dict, ring_list, ring_bond_list)
            bo_sum = np.zeros(atom_num)
            for p, q in bo_dict.keys():
                bo_sum[p] += new_bo_dict[(p, q)]
                bo_sum[q] += new_bo_dict[(p, q)]
            alreadyOctet = (group_list + bo_sum - new_chg_list == 8) & (period_list == 2)
            print("alreadyOctet", np.nonzero(alreadyOctet))
            new_chg_list, new_bo_dict = resolve_chg(
                atom_num,
                bond_num,
                group_list,
                bond_list,
                bond_mapping,
                mve_list,
                neighbor_list,
                chg_mol,
                eIsEven,
                alreadyOctet,
                resolve_step,
            )
            resolve_step += 1

    else:
        new_chg_list, new_bo_dict = chg_list, bo_dict

    bo_matrix2 = np.zeros((atom_num, atom_num))
    for p, q in new_bo_dict.keys():
        bo_matrix2[p][q] = new_bo_dict[(p, q)]
        bo_matrix2[q][p] = new_bo_dict[(p, q)]

    return chg_list, bo_matrix, new_chg_list, bo_matrix2


def compute_chg_and_bo(molecule, chg_mol, resolve=True, cleanUp=True, **kwargs):
    (
        period_list,
        group_list,
        z_list,
        ve_list,
        adj_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        ring_list,
        ring_bond_list,
    ) = get_lists(molecule)

    atom_num, bond_num = len(z_list), len(bond_list)
    eIsEven = int(np.sum(z_list) - chg_mol) % 2 == 0
    resolve_step = 0

    chg_list, bo_dict = maximize_bo(
        atom_num,
        bond_num,
        group_list,
        bond_list,
        bond_mapping,
        ve_list,
        neighbor_list,
        chg_mol,
        eIsEven,
        **kwargs
    )
    # early stop
    if len(bo_dict) == 0:
        return None, None

    # check charge separation
    chg_sep = np.any(chg_list > 0) and np.any(chg_list < 0)

    # charge resolution
    if resolve and chg_sep:
        bo_sum = np.zeros(atom_num)
        for p, q in bo_dict.keys():
            bo_sum[p] += bo_dict[(p, q)]
            bo_sum[q] += bo_dict[(p, q)]
        alreadyOctet = (group_list + bo_sum - chg_list == 8) & (period_list == 2)
        
        # ve_list is now expanded Valence list
        eve_list = get_expanded_list(period_list, ve_list, chg_list, ring_list)
        chg_list, bo_dict = resolve_chg(
            atom_num,
            bond_num,
            group_list,
            bond_list,
            bond_mapping,
            eve_list,
            neighbor_list,
            chg_mol,
            eIsEven,
            alreadyOctet,
            resolve_step
        )

        resolve_step += 1
        
        if cleanUp:
            mve_list = get_modified_list(period_list, eve_list, chg_list, bo_dict, ring_list, ring_bond_list)
            bo_sum = np.zeros(atom_num)
            for p, q in bo_dict.keys():
                bo_sum[p] += bo_dict[(p, q)]
                bo_sum[q] += bo_dict[(p, q)]
            alreadyOctet = (group_list + bo_sum - chg_list == 8) & (period_list == 2)
            chg_list, bo_dict = resolve_chg(
                atom_num,
                bond_num,
                group_list,
                bond_list,
                bond_mapping,
                mve_list,
                neighbor_list,
                chg_mol,
                eIsEven,
                alreadyOctet,
                resolve_step,
            )
            resolve_step += 1

        # error handling
        if len(bo_dict) == 0:
            return None, None

    bo_matrix = np.zeros((atom_num, atom_num))
    for p, q in bo_dict.keys():
        bo_matrix[p][q] = bo_dict[(p, q)]
        bo_matrix[q][p] = bo_dict[(p, q)]

    return chg_list, bo_matrix

