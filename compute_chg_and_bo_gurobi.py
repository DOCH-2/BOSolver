EN_TABLE={
"H":2.300,
"Li":0.912,
"Be":1.576,
"B":2.051,
"C":2.544,
"N":3.066,
"O":3.61,
"F":4.193,
"Ne":4.787,
"Na":0.869,
"Mg":1.293,
"A1":1.613,
"Si":1.916,
"P":2.253,
"S":2.589,
"Cl":2.869,
"Ar":3.242,
"K":0.734,
"Ca":1.034,
"Ga":1.756,
"Ge":1.994,
"As":2.211,
"Se":2.424,
"Br":2.685,
"Kr":2.966,
"Rb":0.706,
"Sr":0.963,
"In":1.656,
"Sn":1.824,
"Sb":1.984,
"Te":2.158,
"I":2.359,
"Xe":2.582
} # J. Am. Chem. Soc. 1989, 111, 25, 9003–9014

from typing import List, Dict, Tuple
from itertools import combinations as comb

from rdkit import Chem

from gurobipy import GRB, Model, LinExpr, Env
from acerxn import chem
import numpy as np


def get_ring_info(z_list, adj_matrix):
    chg_list = np.zeros(len(z_list))
    new_z_list = adj_matrix.sum(axis=-1)
    new_z_list[new_z_list == 1] = 1  # H
    new_z_list[new_z_list == 2] = 8  # O
    new_z_list[new_z_list == 3] = 7  # N
    new_z_list[new_z_list == 4] = 6  # C
    new_z_list[new_z_list == 5] = 15  # P
    new_z_list[new_z_list == 6] = 16  # S
    new_rd = chem.Molecule([new_z_list, adj_matrix, None, chg_list]).get_rd_mol()
    Chem.SanitizeMol(new_rd)

    sssrs = Chem.GetSymmSSSR(new_rd)
    RingInfo = new_rd.GetRingInfo()
    atoms_in_ring = RingInfo.AtomRings()
    bond_rings = RingInfo.BondRings()

    bonds_in_ring = [[] for _ in range(len(sssrs))]
    for ringN, bonds in enumerate(bond_rings):
        for bond in bonds:
            bObj = new_rd.GetBondWithIdx(bond)
            bonds_in_ring[ringN].append((bObj.GetBeginAtomIdx(), bObj.GetEndAtomIdx()))
    #print("bonds in ring", bonds_in_ring)
    #print("bond rings", bond_rings, type(bond_rings))

    ring_neighbors_info = {}

    for aID in set([xx for x in atoms_in_ring for xx in x]):
        atom = new_rd.GetAtomWithIdx(aID)
        ringIDs = RingInfo.AtomMembers(aID)
        ring_bonds= [bond for bond in atom.GetBonds() if bond.IsInRing()]
        ring_dict = dict([(i, []) for i in ringIDs])
        for bond in ring_bonds:
            for bRID in RingInfo.BondMembers(bond.GetIdx()):
                ring_dict[bRID].append((bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()))
        ring_neighbors_info[aID] = ring_dict.values()

    """
    # Example: Benzene
    >>> atoms_in_ring: [[0, 1, 2, 3, 4, 5]]
    >>> bonds_in_ring: [[(0,1),(1,2),(2,3),(3,4),(4,5),(5,0)]]
    >>> ring_neighbors_info: {0: [[(0, 5), (0, 1)]], 
                              1: [[(1, 2), (0, 1)]], 
                              2: [[(1, 2), (2, 3)]], 
                              3: [[(3, 4), (2, 3)]], 
                              4: [[(3, 4), (4, 5)]], 
                              5: [[(0, 5), (4, 5)]]}
    """

    return atoms_in_ring, bonds_in_ring, ring_neighbors_info


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
        elif period_list[i] == 2:
            # not considering expanded octet here
            ve_list[i] = 8
        else:
            ve_list[i] = 18

    # ring membership
    _, _, ring_neighbors_info = get_ring_info(z_list, adj_matrix)

    # electronegativity
    en_list = np.array([EN_TABLE[chem.periodic_table[z-1]] for z in z_list])

    return (
        period_list,
        group_list,
        z_list,
        ve_list,
        adj_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        ring_neighbors_info,
        en_list
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


def get_modified_list(
    period_list, eve_list, chg_list, bo_dict, ring_list, ring_bond_list
):
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
    en_list,
    ring_neighbors_info,
    chg_mol,
    eIsEven,
    **kwargs,
):
    # early stop
    if atom_num == 1:
        return np.array([chg_mol]), {}

    ### model construction
    env = Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.start()
    model = Model("maximize_bo", env=env)
    verbose = kwargs.get("printOptLog", False)
    Xsingle = kwargs.get("HalogenConstraint", False)
    cleanUp = kwargs.get("cleanUp", False) and (len(ring_neighbors_info) > 0)

    # Bond order
    # db: double bond flag
    # tb: triple bond flag
    # the bond order would be 1, 2, 3
    # from the nature, the constraint db + tb <= 1 should be given
    # and bond order is represented as (1 + db + 2 * tb)
    db = model.addVars(bond_num, name="dbFlag", vtype=GRB.BINARY)
    tb = model.addVars(bond_num, name="tbFlag", vtype=GRB.BINARY)

    model.addConstrs(
        (db[i] + tb[i] <= 1 for i in range(bond_num)), name="BondOrderFlag"
    )

    # t1: formal charge
    # t1[2i]: fc+ | t1[2i+1]: fc-
    # t1[2i] - t1[2i+1] : fc of atom i
    # t1[2i] + t1[2i+1] : abs(fc) of atom i
    t1 = model.addVars(2 * atom_num, lb=0, name="t1", vtype=GRB.INTEGER)

    # t2: formal charge for weighted objective function
    # weight considering electronegativity
    # t2 = model.addVars(2 * atom_num, name="t2", vtype=GRB.CONTINUOUS)

    # Halogen Constraints
    # Halogen atoms, especially F and Cl, are not allowed to have
    # formal charge whose absolute value is greater than 0,
    # and they have only one single bond.
    # Br occassionaly has +1 formal charge, with 2 single bonds (as a three-membered ring)
    # But for now, just simply apply the rule to all halogen atoms.
    if Xsingle:
        for i in np.where(group_list == 7)[0]:
            model.addConstr(t1[2 * i] == 0, name=f"halogen+_{i}")
            model.addConstr(t1[2 * i + 1] == 0, name=f"halogen-_{i}")
    #   model.addConstr(t1[2 * i] + t1[2 * i + 1] == 0, name=f"halogen_{i}")

    # Ring Constraints
    # the atom in the ring should have at most one double bond.
    # if cleanUp:
    #    ring_members = np.unique(sum(ring_list, [])).tolist()
    #    ring_constrs = [
    #        LinExpr() for _ in range(len(ring_members))
    #    ]  # same order with ring_members
    #    # print("rings", ring_bond_list)
    #    for ring_bond in ring_bond_list:
    #        for i, j in ring_bond:
    #            # ring_memb_idx1 =
    #            # ring_memb_idx2 =
    #            ring_constrs[ring_members.index(i)].add(bo[bond_mapping[(i, j)]])
    #            ring_constrs[ring_members.index(j)].add(bo[bond_mapping[(i, j)]])
    #    for i, ring_constr in enumerate(ring_constrs):
    #        model.addConstr(
    #            ring_constr <= ring_constr.size() + 1,
    #            name=f"ring_{ring_members[i]}",
    #        )

    # even: dummy variable to force no. of electrons even
    even = model.addVars(atom_num, name="even", vtype=GRB.INTEGER)

    ### constraints construction
    chg_constr = LinExpr()  # charge conservation rule
    for i in range(atom_num):
        lp_constr = LinExpr()  # lone pair rule
        ve_constr = LinExpr()  # valence rule
        X_constr = (
            LinExpr() if (Xsingle and group_list[i] == 7) else None
        )  # Halogen Constraint

        chg_constr.add(t1[2 * i] - t1[2 * i + 1])
        lp_constr.add(t1[2 * i] - t1[2 * i + 1])
        ve_constr.add(-t1[2 * i] + t1[2 * i + 1])

        # summation over bond
        for j in neighbor_list[i]:
            a, b = i, j
            if a > b:
                a, b = b, a

            bo = LinExpr(
                1 + db[bond_mapping[(a, b)]] + tb[bond_mapping[(a, b)]] * 2
            )  # bond order

            lp_constr.add(bo)
            ve_constr.add(bo)
            if X_constr is not None:
                X_constr.add(bo)

        # the number of lone pair should not be negative
        model.addConstr(lp_constr <= group_list[i], name=f"lp_{i}")

        # the number of valence electron is not greater than maximum valence
        model.addConstr(ve_constr + group_list[i] <= ve_list[i], name=f"ve_{i}")

        if eIsEven:
            # the number of valence electron is even number (no radical rule!)
            model.addConstr(ve_constr + group_list[i] == 2 * even[i], name=f"noRad_{i}")

        # halogen atoms have only one single bond
        if X_constr is not None:
            model.addConstr(X_constr == 1, name=f"X_{i}")

        ### TODO: Ring Constraint
        if cleanUp and (i in ring_neighbors_info):
            for n, neighbor_bond_list in enumerate(ring_neighbors_info[i]):
                ring_constr = LinExpr()  # ring constraint
                for ring_bond in neighbor_bond_list:
                    ring_constr.add(
                        db[bond_mapping[ring_bond]] + tb[bond_mapping[ring_bond]]
                    )
                model.addConstr(ring_constr <= 1, name=f"ring_{i}_{n}")

    model.addConstr(chg_constr == chg_mol, name="chg_consv")

    ### optimization
    max_bo_obj = LinExpr()  # bond maximization
    min_fc_obj = LinExpr()  # formal charge minimization
    min_en_obj = LinExpr()  # electronegativity minimization

    # bond maximization
    for i in range(bond_num):
        bo = LinExpr(1 + db[i] + tb[i] * 2)  # bond order
        max_bo_obj.add(bo)
    # formal charge minimization
    for i in range(atom_num):
        min_fc_obj.add(t1[2 * i] + t1[2 * i + 1])
        # min_wfc_obj.add((0.1 * group_list[i] + 0.6) * (t2[2 * i] + t2[2 * i + 1]))
        
    for i, en in enumerate(en_list):
        min_en_obj.add(en * (t1[2 * i] - t1[2 * i + 1]))

    bo_priority = 3 # bond order maximization priority
    chg_priority = 2 # charge separation priority
    en_priority = 1 # electronegativity priority


    if kwargs.get("mode", "") == "fc":
        bo_priority, chg_priority = chg_priority, bo_priority

    model.setObjectiveN(
        max_bo_obj, 1, priority=bo_priority, weight=GRB.MAXIMIZE, name="max_bo"
    )
    model.setObjectiveN(
        min_fc_obj, 2, priority=chg_priority, weight=GRB.MINIMIZE, name="min_fc"
    )
    model.setObjectiveN(
        min_en_obj, 3, priority=en_priority, weight=GRB.MINIMIZE, name="min_en"
    )


    # Gurobi optimization
    # model.setParam(GRB.Param.OutputFlag, 0)
    model.setParam(GRB.Param.TimeLimit, 1)
    if verbose:
        model.write("record.lp")

    model.optimize()

    # error handling
    if model.status != GRB.Status.OPTIMAL:
        # model.write("output.bas")
        # model.write("output.mst")
        # model.write("output.sol")
        return np.zeros(atom_num), {}

    # result record
    if verbose:
        nSol = model.SolCount
        for s in range(nSol):
            model.params.SolutionNumber = s
            model.write(f"output{s}.sol")

    # retrieval
    bo_dict = {}
    chg_list = np.zeros(atom_num, dtype=np.int64)
    for i in range(bond_num):
        bo = 1 + int(db[i].X) + 2 * int(tb[i].X)
        bo_dict[bond_list[i]] = int(bo)
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
    **kwargs,
):
    if atom_num == 1:
        return np.array([chg_mol]), {}

    ### model construction
    env = Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.start()
    verbose = kwargs.get("printOptLog", False)
    Xsingle = kwargs.get("HalogenConstraint", False)

    model = Model(f"resolve_chg{stepIdx}", env=env)

    # bo: bond order
    bo = model.addVars(bond_num, lb=1, name="bo", vtype=GRB.INTEGER)

    # t1: formal charge
    t1 = model.addVars(2 * atom_num, lb=0, name="t1", vtype=GRB.INTEGER)
    # t2: formal charge for weighted objective function
    # weight considering electronegativity
    # t2 = model.addVars(2 * atom_num, name="t2", vtype=GRB.CONTINUOUS)

    # Halogen atoms, especially F and Cl, are not allowed to have
    # formal charge whose absolute value is greater than 0,
    # and they have only one single bond.
    # Br occassionaly has +1 formal charge, with 2 single bonds (as a three-membered ring)
    # But for now, just simply apply the rule to all halogen atoms.
    if Xsingle:
        for i in np.where(group_list == 7)[0]:
            model.addConstr(t1[2 * i] == 0, name=f"halogen+_{i}")
            model.addConstr(t1[2 * i + 1] == 0, name=f"halogen-_{i}")
            # model.addConstr(t1[2 * i] + t1[2 * i + 1] == 0, name=f"halogen_{i}")

    # even: dummy variable to force no. of electrons even
    even = model.addVars(atom_num, name="even", vtype=GRB.INTEGER)

    ### constraints construction
    chg_constr = LinExpr()  # charge conservation rule
    octet_constr = LinExpr()  #

    for i in range(atom_num):
        lp_constr = LinExpr()  # lone pair rule
        eve_constr = LinExpr()  # expanded valence rule
        X_constr = (
            LinExpr() if (Xsingle and group_list[i] == 7) else None
        )  # Halogen rule

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
            if X_constr is not None:
                X_constr.add(bo[bond_mapping[(a, b)]])

        model.addConstr(lp_constr <= group_list[i], name=f"lp_{i}")
        if bool(alreadyOctet[i]):
            model.addConstr(eve_constr + group_list[i] == 8, name=f"eve_{i}")
        else:
            model.addConstr(eve_constr <= eve_list[i] - group_list[i], name=f"eve_{i}")

        # halogen atoms have only one single bond
        if X_constr is not None:
            model.addConstr(X_constr == 1, name=f"X_{i}")

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
    # model.setParam(GRB.Param.OutputFlag, 0)
    model.setParam(GRB.Param.TimeLimit, 1)
    if verbose:
        model.write(f"record_chg{stepIdx}.lp")

    model.optimize()

    # error handling
    if model.status != GRB.Status.OPTIMAL:
        return np.zeros(atom_num), {}

    # result record
    nSol = model.SolCount
    # print("nSol", nSol)
    for s in range(nSol):
        model.params.SolutionNumber = s
        if verbose:
            model.write(f"output_chg{stepIdx}_{s}.sol")

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

    # def compute_chg_and_bo_debug(molecule, chg_mol, resolve=True, cleanUp=True, **kwargs):
    #    (
    #        period_list,
    #        group_list,
    #        z_list,
    #        ve_list,
    #        adj_list,
    #        bond_list,
    #        bond_mapping,
    #        neighbor_list,
    #        ring_neighbors_info,
    #    ) = get_lists(molecule)
    #    atom_num, bond_num = len(z_list), len(bond_list)
    #    eIsEven = int(np.sum(z_list) - chg_mol) % 2 == 0
    #    resolve_step = 0
    #
    #    chg_list, bo_dict = maximize_bo(
    #        atom_num,
    #        bond_num,
    #        group_list,
    #        bond_list,
    #        bond_mapping,
    #        ve_list,
    #        neighbor_list,
    #        chg_mol,
    #        eIsEven,
    #        **kwargs,
    #    )
    #    # early stop
    #    if len(bo_dict) == 0:
    #        return None, None, None, None
    #
    #    bo_matrix = np.zeros((atom_num, atom_num))
    #    for p, q in bo_dict.keys():
    #        bo_matrix[p][q] = bo_dict[(p, q)]
    #        bo_matrix[q][p] = bo_dict[(p, q)]
    #
    #    # check charge separation
    #    chg_sep = np.any(chg_list > 0) and np.any(chg_list < 0)
    #
    #    # charge resolution
    #    if resolve and chg_sep:
    #        bo_sum = np.zeros(atom_num)
    #        for p, q in bo_dict.keys():
    #            bo_sum[p] += bo_dict[(p, q)]
    #            bo_sum[q] += bo_dict[(p, q)]
    #        alreadyOctet = (group_list + bo_sum - chg_list == 8) & (period_list == 2)
    #        print("alreadyOctet", np.nonzero(alreadyOctet))
    #
    #        # ve_list is now expanded Valence list
    #        eve_list = get_expanded_list(period_list, ve_list, chg_list, ring_list)
    #        new_chg_list, new_bo_dict = resolve_chg(
    #            atom_num,
    #            bond_num,
    #            group_list,
    #            bond_list,
    #            bond_mapping,
    #            eve_list,
    #            neighbor_list,
    #            chg_mol,
    #            eIsEven,
    #            alreadyOctet,
    #            resolve_step,
    #        )
    #        resolve_step += 1
    #
    #        if cleanUp:
    #            mve_list = get_modified_list(
    #                period_list,
    #                eve_list,
    #                new_chg_list,
    #                new_bo_dict,
    #                ring_list,
    #                ring_bond_list,
    #            )
    #            bo_sum = np.zeros(atom_num)
    #            for p, q in bo_dict.keys():
    #                bo_sum[p] += new_bo_dict[(p, q)]
    #                bo_sum[q] += new_bo_dict[(p, q)]
    #            alreadyOctet = (group_list + bo_sum - new_chg_list == 8) & (
    #                period_list == 2
    #            )
    #            print("alreadyOctet", np.nonzero(alreadyOctet))
    #            new_chg_list, new_bo_dict = resolve_chg(
    #                atom_num,
    #                bond_num,
    #                group_list,
    #                bond_list,
    #                bond_mapping,
    #                mve_list,
    #                neighbor_list,
    #                chg_mol,
    #                eIsEven,
    #                alreadyOctet,
    #                resolve_step,
    #            )
    #            resolve_step += 1
    #
    #    else:
    #        new_chg_list, new_bo_dict = chg_list, bo_dict
    #
    #    bo_matrix2 = np.zeros((atom_num, atom_num))
    #    for p, q in new_bo_dict.keys():
    #        bo_matrix2[p][q] = new_bo_dict[(p, q)]
    #        bo_matrix2[q][p] = new_bo_dict[(p, q)]
    #
    #    return chg_list, bo_matrix, new_chg_list, bo_matrix2


def compute_chg_and_bo(molecule, chg_mol, resolve=False, cleanUp=True, **kwargs):
    (
        period_list,
        group_list,
        z_list,
        ve_list,
        adj_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        ring_neighbors_info,
        en_list,
    ) = get_lists(molecule)

    atom_num, bond_num = len(z_list), len(bond_list)
    eIsEven = int(np.sum(z_list) - chg_mol) % 2 == 0
    resolve_step = 0
    kwargs["cleanUp"] = cleanUp

    chg_list, bo_dict = maximize_bo(
        atom_num,
        bond_num,
        group_list,
        bond_list,
        bond_mapping,
        ve_list,
        neighbor_list,
        en_list,
        ring_neighbors_info,
        chg_mol,
        eIsEven,
        **kwargs,
    )

    # early stop
    if len(bo_dict) == 0:
        return None, None

    # check charge separation
    chg_sep = np.any(chg_list > 0) and np.any(chg_list < 0)

    # charge resolution
    if resolve and chg_sep:
        print("You shall not pass here")
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
            resolve_step,
            **kwargs,
        )

        resolve_step += 1

        # error handling
        # if len(bo_dict) == 0:
        #    return None, None

        if cleanUp:
            mve_list = get_modified_list(
                period_list, ve_list, chg_list, bo_dict, ring_list, ring_bond_list
            )
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
                **kwargs,
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
