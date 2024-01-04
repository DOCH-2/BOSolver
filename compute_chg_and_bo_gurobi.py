from typing import List, Dict, Tuple

from gurobipy import GRB, Model, LinExpr
from acerxn import chem
import numpy as np


def get_lists(molecule: chem.Molecule):
    period_list, group_list = molecule.get_period_group_list()
    adj_matrix = np.copy(molecule.get_matrix("adj"))
    adj_list = np.sum(adj_matrix, axis=1)

    new_adj = np.copy(adj_matrix)
    new_adj[np.tril_indices_from(new_adj)] = 0
    neighbor_index = np.nonzero(new_adj)

    reduced_neighbor_dict = {}
    for i, k in enumerate(neighbor_index[0]):
        if k not in reduced_neighbor_dict:
            reduced_neighbor_dict[k] = [neighbor_index[1][i]]
        else:
            reduced_neighbor_dict[k].append(neighbor_index[1][i])

    ve_list = np.zeros_like(group_list)
    z_list = molecule.get_z_list()
    for i in range(len(group_list)):
        if period_list[i] == 1:
            ve_list[i] = 2
        else:
            # not considering expanded octet here
            ve_list[i] = 8

    bond_list = molecule.get_bond_list(False)
    bond_mapping = {key: val for val, key in enumerate(bond_list)}
    neighbor_list = molecule.get_neighbor_list()
    return (
        period_list,
        group_list,
        z_list,
        ve_list,
        adj_list,
        bond_list,
        bond_mapping,
        neighbor_list,
    )


def get_extended_lists(period_list, ve_list, chg_list):
    extended_idx = period_list > 2

    eve_list = np.copy(ve_list)
    eve_list[extended_idx] += 2*np.where(chg_list > 0, chg_list, 0)[extended_idx]
    print("ve_list", ve_list, len(ve_list))
    print("eve_list", eve_list, len(eve_list))

    return eve_list


def maximize_bo(
    atom_num,
    bond_num,
    group_list,
    bond_list,
    bond_mapping,
    ve_list,
    neighbor_list,
    chg_mol,
    mode,
):
    # early stop
    if atom_num == 1:
        return np.array([chg_mol]), {}

    ### model construction
    model = Model("maximize_bo")

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
        model.addConstr(ve_constr <= ve_list[i] - group_list[i], name=f"ve_{i}")
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

    model.setObjectiveN(max_bo_obj, 0, priority=2, weight=GRB.MAXIMIZE, name="max_bo")
    model.setObjectiveN(min_fc_obj, 1, priority=1, weight=GRB.MINIMIZE, name="min_fc")
    # if mode == "heuristics":
    #    model.setObjectiveN(min_wfc_obj, 2, 0, GRB.MINIMIZE)

    # Gurobi optimization
    # model.setParam(GRB.Param.OutputFlag, 0)
    model.setParam(GRB.Param.TimeLimit, 1)
    model.write("record.lp")

    model.optimize()

    # error handling
    if model.status != GRB.Status.OPTIMAL:
        return np.zeros(atom_num), {}

    # result record
    nSol = model.SolCount
    print("nSol", nSol)
    for s in range(nSol):
        model.params.SolutionNumber = s
        model.write(f"output{s}.sol")

    # retrieval
    bo_dict = {}
    chg_list = np.zeros(atom_num, dtype=np.int64)
    for i in range(bond_num):
        bo_dict[bond_list[i]] = int(bo[i].X)
    for i in range(atom_num):
        chg_list[i] = int(t1[2 * i].X) - int(t1[2 * i + 1].X)

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
    mode,
):
    if atom_num == 1:
        return np.array([chg_mol]), {}

    model = Model("resolve_chg")

    # bo: bond order
    bo = model.addVars(bond_num, lb=1, name="bo", vtype=GRB.INTEGER)

    # t1: formal charge
    t1 = model.addVars(2 * atom_num, lb=0, name="t1", vtype=GRB.INTEGER)
    # t2: formal charge for weighted objective function
    # weight considering electronegativity
    # t2 = model.addVars(2 * atom_num, name="t2", vtype=GRB.CONTINUOUS)

    ### constraints construction
    chg_constr = LinExpr()  # charge conservation rule
    for i in range(atom_num):
        lp_constr = LinExpr()  # lone pair rule
        eve_constr = LinExpr()  # extended valence rule

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
        model.addConstr(eve_constr <= eve_list[i] - group_list[i], name=f"eve_{i}")
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
    # model.setParam(GRB.Param.OutputFlag, 0)
    model.setParam(GRB.Param.TimeLimit, 1)
    model.write("record_chg.lp")

    model.optimize()

    # error handling
    if model.status != GRB.Status.OPTIMAL:
        return np.zeros(atom_num), {}

    # result record
    nSol = model.SolCount
    print("nSol", nSol)
    for s in range(nSol):
        model.params.SolutionNumber = s
        model.write(f"output_chg{s}.sol")

    # retrieval
    bo_dict = {}
    chg_list = np.zeros(atom_num)
    for i in range(bond_num):
        bo_dict[bond_list[i]] = int(bo[i].X)
    for i in range(atom_num):
        chg_list[i] = int(t1[2 * i].X) - int(t1[2 * i + 1].X)
    return chg_list, bo_dict


def compute_chg_and_bo(molecule, chg_mol, resolve=True, mode="heuristics"):
    (
        period_list,
        group_list,
        z_list,
        ve_list,
        adj_list,
        bond_list,
        bond_mapping,
        neighbor_list,
    ) = get_lists(molecule)
    atom_num, bond_num = len(z_list), len(bond_list)

    chg_list, bo_dict = maximize_bo(
        atom_num,
        bond_num,
        group_list,
        bond_list,
        bond_mapping,
        ve_list,
        neighbor_list,
        chg_mol,
        mode,
    )
    # early stop
    if len(bo_dict) == 0:
        return None, None, None, None

    bo_matrix = np.zeros((atom_num, atom_num))
    for p, q in bo_dict.keys():
        bo_matrix[p][q] = bo_dict[(p, q)]
        bo_matrix[q][p] = bo_dict[(p, q)]

    # charge resolution
    if resolve:
        eve_list = get_extended_lists(period_list, ve_list, chg_list)
        new_chg_list, new_bo_dict = resolve_chg(
            atom_num,
            bond_num,
            group_list,
            bond_list,
            bond_mapping,
            eve_list,
            neighbor_list,
            chg_mol,
            mode,
        )
    else:
        return chg_list, bo_matrix, None, None

    bo_matrix2 = np.zeros((atom_num, atom_num))
    for p, q in new_bo_dict.keys():
        bo_matrix2[p][q] = new_bo_dict[(p, q)]
        bo_matrix2[q][p] = new_bo_dict[(p, q)]

    return chg_list, bo_matrix, new_chg_list, bo_matrix2
