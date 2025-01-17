from typing import List, Dict, Tuple
from itertools import combinations as comb
import sys

from rdkit import Chem
from sqlalchemy import over

from acerxn import chem
from gurobipy import GRB, Model, LinExpr, Env
import numpy as np
from scipy import spatial

EN_TABLE = {
    "H": 2.300,
    "Li": 0.912,
    "Be": 1.576,
    "B": 2.051,
    "C": 2.544,
    "N": 3.066,
    "O": 3.61,
    "F": 4.193,
    "Ne": 4.787,
    "Na": 0.869,
    "Mg": 1.293,
    "A1": 1.613,
    "Si": 1.916,
    "P": 2.253,
    "S": 2.589,
    "Cl": 2.869,
    "Ar": 3.242,
    "K": 0.734,
    "Ca": 1.034,
    "Ga": 1.756,
    "Ge": 1.994,
    "As": 2.211,
    "Se": 2.424,
    "Br": 2.685,
    "Kr": 2.966,
    "Rb": 0.706,
    "Sr": 0.963,
    "In": 1.656,
    "Sn": 1.824,
    "Sb": 1.984,
    "Te": 2.158,
    "I": 2.359,
    "Xe": 2.582,
}  # J. Am. Chem. Soc. 1989, 111, 25, 9003–9014


def get_gurobi_model_env():
    env = Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.start()
    model = Model("compute_chg_and_bo", env=env)

    return model, env


def get_adj_matrix_from_distance3(molecule):
    n = len(molecule.atom_list)
    radius_list = molecule.get_radius_list()
    radius_matrix = np.repeat(radius_list, n).reshape((n, n))
    criteria_matrix = (radius_matrix + radius_matrix.T) + np.ones(
        (n, n)
    ) * 0.45  # J. Chem. Inf. Comput. Sci. 1992, 32, 401-406
    coordinate_list = molecule.get_coordinate_list()
    distance_matrix = spatial.distance_matrix(coordinate_list, coordinate_list)
    adj = np.where(distance_matrix < criteria_matrix, 1, 0)
    np.fill_diagonal(adj, 0)
    return adj


def get_adj_matrix_from_distance4(molecule, coeff=1.15):
    n = len(molecule.atom_list)
    radius_list = molecule.get_radius_list()
    radius_matrix = np.repeat(radius_list, n).reshape((n, n))
    criteria_matrix = (radius_matrix + radius_matrix.T) * coeff
    coordinate_list = molecule.get_coordinate_list()
    distance_matrix = spatial.distance_matrix(coordinate_list, coordinate_list)
    adj = np.where(
        ((distance_matrix < criteria_matrix) | (distance_matrix < 0.80)), 1, 0
    )
    np.fill_diagonal(adj, 0)
    return adj


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
    # print("bonds in ring", bonds_in_ring)
    # print("bond rings", bond_rings, type(bond_rings))

    ring_neighbors_info = {}

    for aID in set([xx for x in atoms_in_ring for xx in x]):
        atom = new_rd.GetAtomWithIdx(aID)
        ringIDs = RingInfo.AtomMembers(aID)
        ring_bonds = [bond for bond in atom.GetBonds() if bond.IsInRing()]
        ring_dict = dict([(i, []) for i in ringIDs])
        for bond in ring_bonds:
            for bRID in RingInfo.BondMembers(bond.GetIdx()):
                ring_dict[bRID].append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
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


def get_lists(molecule: chem.Molecule, strict=True):
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
            ve_list[i] = 8 if strict else 2 * group_list[i]
            # ve_list[i] = 18

    # ring membership
    _, _, ring_neighbors_info = get_ring_info(z_list, adj_matrix)

    # electronegativity
    en_list = np.array([EN_TABLE[chem.periodic_table[z - 1]] for z in z_list])

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
        en_list,
    )


def get_expanded_ve_list(period_list, group_list, ve_list, chg_list):
    # apply expanded octet rule only if satisfies the following conditions
    # 1. period > 2
    # 2. has non-zero formal charge
    # Ring constraint is not considered here anymore

    # ring_members = np.unique(sum(ring_list, []))
    # in_ring = np.zeros_like(period_list)
    # if len(ring_members) > 0:
    #    in_ring[ring_members] = 1
    # in_ring = in_ring.astype(bool)

    # expanded_idx = (period_list > 2) & ~in_ring
    expanded_idx = period_list > 2
    eve_list = np.copy(ve_list)
    eve_list[expanded_idx] += (
        2 * np.where(chg_list > 0, np.minimum(group_list, chg_list), 0)[expanded_idx]
    )

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


# TODO: Code Acceleration
# Use single Gurobi model for both optimization and resolution


def maximize_bo(
    atom_num,
    bond_num,
    period_list,
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
        return np.array([chg_mol]), {}, (None, None, None)

    ### model construction
    env = Env(empty=True)
    env.setParam("OutputFlag", 0)
    # env.setParam("DualReductions", 0)
    # env.setParam("LogFile", "gurobi.log")
    env.start()
    model = Model("maximize_bo", env=env)
    verbose = kwargs.get("printOptLog", False)
    Xsingle = kwargs.get("HalogenConstraint", False)
    cleanUp = kwargs.get("cleanUp", False) and (len(ring_neighbors_info) > 0)
    M_list = kwargs.get("MetalCenters", [])

    db_starts = kwargs.get("db_starts", [0] * bond_num)
    tb_starts = kwargs.get("tb_starts", [0] * bond_num)
    t1_starts = kwargs.get("t1_starts", [0] * 2 * atom_num)

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
    # b1 = model.addVars(bond_num, lb=1, ub=3, name="b", vtype=GRB.INTEGER)

    # t1: formal charge
    # t1[2i]: fc+ | t1[2i+1]: fc-
    # t1[2i] - t1[2i+1] : fc of atom i
    # t1[2i] + t1[2i+1] : abs(fc) of atom i
    t1 = model.addVars(2 * atom_num, lb=0, name="t1", vtype=GRB.INTEGER)

    # t2: formal charge for weighted objective function
    # weight considering electronegativity
    # t2 = model.addVars(2 * atom_num, name="t2", vtype=GRB.CONTINUOUS)

    # o: octet distance
    # the distance between the number of valence elctrons and
    # the octet number(2 for 1st period, 8 for 2nd or higher period)
    # o = 8 - (g - c + b)

    # Set Initial Values
    for i in range(bond_num):
        db[i].Start = db_starts[i]
        tb[i].Start = tb_starts[i]
    # for i in range(bond_num):
    #    b1[i].Start = 1 + db_starts[i] + 2 * tb_starts[i]
    for i in range(2 * atom_num):
        t1[i].Start = t1_starts[i]

    # TODO: Revision of Halogen Constraint
    # Halogen atoms, especially Cl and Br, are not allowed for
    # following the extended octet rule.
    # RDKit does not allow Cl and Br to have valence state greater than 1

    # even: dummy variable to force no. of electrons even
    even = model.addVars(atom_num, name="even", vtype=GRB.INTEGER)

    ### objectives and constraints construction
    # objective functions
    min_od_obj = LinExpr()  # octet distance minimization (octet maximization)
    min_fc_obj = LinExpr()  # formal charge minimization
    max_bo_obj = LinExpr()  # bond maximization
    min_en_obj = LinExpr()  # electronegativity minimization

    # constraints
    chg_constr = LinExpr()  # charge conservation rule

    for i in range(atom_num):
        lp_constr = LinExpr()  # lone pair rule
        ve_constr = LinExpr()  # valence electron
        # X_constr = (
        #    LinExpr()
        #    if (Xsingle and group_list[i] == 7 and period_list[i] <= 4)
        #    else None
        # )  # Halogen Constraint

        ve_constr.addConstant(group_list[i])

        chg_constr.add(t1[2 * i] - t1[2 * i + 1])
        lp_constr.add(t1[2 * i] - t1[2 * i + 1])
        ve_constr.add(-t1[2 * i] + t1[2 * i + 1])
        min_fc_obj.add(t1[2 * i] + t1[2 * i + 1])
        min_en_obj.add(en_list[i] * (t1[2 * i] - t1[2 * i + 1]))

        # summation over bond
        for j in neighbor_list[i]:
            a, b = i, j
            if a > b:
                a, b = b, a

            bo = LinExpr(
                1 + db[bond_mapping[(a, b)]] + tb[bond_mapping[(a, b)]] * 2
            )  # bond order

            # bo = b1[bond_mapping[(a, b)]]
            lp_constr.add(bo)
            ve_constr.add(bo)
            # if X_constr is not None:
            #    X_constr.add(bo)

            max_bo_obj.add(bo)

            if Xsingle and group_list[i] == 7 and period_list[i] <= 4:
                model.addConstr(bo == 1, name=f"XC_{i}")

            if i in M_list:
                model.addConstr(bo == 1, name=f"SB_{i}_{j}")

        # the number of lone pair should not be negative
        model.addConstr(lp_constr <= group_list[i], name=f"lp_{i}")

        # TODO: octet rule
        # octet distance
        if period_list[i] == 1:
            min_od_obj.add(2 - ve_constr)
            model.addConstr(2 - ve_constr >= 0, name=f"od_{i}")
        elif period_list[i] == 2 or len(neighbor_list[i]) <= 4:
            min_od_obj.add(8 - ve_constr)
            model.addConstr(8 - ve_constr >= 0, name=f"od_{i}")

        # the number of valence electron is not greater than maximum valence
        # model.addConstr(ve_constr + group_list[i] <= ve_list[i], name=f"ve_{i}")
        #
        # the number of valence electron is equal to maximum valence (8 for 2nd period atoms)
        # model.addConstr(ve_constr + group_list[i] == ve_list[i], name=f"ve_{i}")

        if eIsEven:
            # the number of valence electron is even number (no radical rule!)
            model.addConstr(ve_constr == 2 * even[i], name=f"noRad_{i}")

        # halogen atoms have only one single bond
        # halogens might have any bond (halogen anions), and in such case, does not apply the constraint
        # if X_constr is not None and X_constr.size() > 0:
        #    model.addConstr(X_constr == 1, name=f"X_{i}")

        # Ring Constraint
        if cleanUp and (i in ring_neighbors_info):
            for n, neighbor_bond_list in enumerate(ring_neighbors_info[i]):
                ring_constr = LinExpr()  # ring constraint
                for ring_bond in neighbor_bond_list:
                    ring_constr.add(
                        db[bond_mapping[ring_bond]] + tb[bond_mapping[ring_bond]]
                    )
                model.addConstr(ring_constr <= 1, name=f"ring_{i}_{n}")

    model.addConstr(chg_constr == chg_mol, name="chg_consv")

    ## formal charge minimization
    # for i in range(atom_num):
    #    min_fc_obj.add(t1[2 * i] + t1[2 * i + 1])
    #    # min_wfc_obj.add((0.1 * group_list[i] + 0.6) * (t2[2 * i] + t2[2 * i + 1]))
    ## bond maximization
    # for i in range(bond_num):
    #    bo = LinExpr(1 + db[i] + tb[i] * 2)  # bond order
    #    max_bo_obj.add(bo)
    ## electronegativity minimization
    # for i, en in enumerate(en_list):
    #    min_en_obj.add(en * (t1[2 * i] - t1[2 * i + 1]))

    od_priority = 4  # octet distance priority
    bo_priority = 3  # bond order maximization priority
    chg_priority = 2  # charge separation priority
    en_priority = 1  # electronegativity priority

    if kwargs.get("mode", "") == "fc":
        bo_priority, chg_priority = chg_priority, bo_priority

    model.setObjectiveN(
        min_od_obj, 1, priority=od_priority, weight=GRB.MINIMIZE, name="min_od"
    )

    model.setObjectiveN(
        max_bo_obj, 2, priority=bo_priority, weight=GRB.MAXIMIZE, name="max_bo"
    )
    model.setObjectiveN(
        min_fc_obj, 3, priority=chg_priority, weight=GRB.MINIMIZE, name="min_fc"
    )
    # model.setObjectiveN(
    #   min_en_obj, 4, priority=en_priority, weight=GRB.MINIMIZE, name="min_en"
    # )

    # Gurobi optimization
    # model.setParam(GRB.Param.OutputFlag, 0)
    # model.setParam(GRB.Param.TimeLimit, 1)
    if verbose:
        model.write("record.lp")
    # model.write("model.lp")
    # model.computeIIS()
    # model.write("model.ilp")
    model.optimize()

    # error handling
    if model.status != GRB.Status.OPTIMAL:
        # model.status=3 means infeasible
        print(
            f"maximize_bo: Optimization failed. (status: {model.status})",
            file=sys.stderr,
        )
        # model.write("record.lp")
        # model.write("output.bas")
        # model.write("output.mst")
        # model.write("output.sol")
        return None, None, (None, None, None)

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

    db_values = [int(db[i].X) for i in range(bond_num)]
    tb_values = [int(tb[i].X) for i in range(bond_num)]
    t1_values = [int(t1[i].X) for i in range(2 * atom_num)]
    model.close()
    env.close()

    return chg_list, bo_dict, (db_values, tb_values, t1_values)


def resolve_chg(
    atom_num,
    bond_num,
    period_list,
    group_list,
    bond_list,
    bond_mapping,
    eve_list,
    neighbor_list,
    en_list,
    ring_neighbors_info,
    chg_mol,
    eIsEven,
    overcharged,
    db_starts,
    tb_starts,
    t1_starts,
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
    cleanUp = kwargs.get("cleanUp", False) and (len(ring_neighbors_info) > 0)
    M_list = kwargs.get("MetalCenters", [])

    model = Model(f"resolve_chg{stepIdx}", env=env)

    # bo: bond order
    db = model.addVars(bond_num, name="dbFlag", vtype=GRB.BINARY)
    tb = model.addVars(bond_num, name="tbFlag", vtype=GRB.BINARY)
    model.addConstrs(
        (db[i] + tb[i] <= 1 for i in range(bond_num)), name="BondOrderFlag"
    )

    # t1: formal charge
    t1 = model.addVars(2 * atom_num, lb=0, name="t1", vtype=GRB.INTEGER)
    # import numpy as npt2: formal charge for weighted objective function
    # weight considering electronegativity
    # t2 = model.addVars(2 * atom_num, name="t2", vtype=GRB.CONTINUOUS)

    # Set Initial Values
    for i in range(bond_num):
        db[i].Start = db_starts[i]
        tb[i].Start = tb_starts[i]
    for i in range(2 * atom_num):
        t1[i].Start = t1_starts[i]

    # TODO: Revision of Halogen Constraint
    # Halogen atoms, especially Cl and Br, are not allowed for
    # following the extended octet rule.
    # RDKit does not allow Cl and Br to have valence state greater than 1

    # even: dummy variable to force no. of electrons even
    even = model.addVars(atom_num, name="even", vtype=GRB.INTEGER)

    ### constraints construction
    chg_constr = LinExpr()  # charge conservation rule

    model.update()

    for i in range(atom_num):
        lp_constr = LinExpr()  # lone pair rule
        ve_constr = LinExpr()  # valence electron
        # X_constr = (
        #    LinExpr()
        #    if (Xsingle and group_list[i] == 7 and period_list[i] <= 4)
        #    else None
        # )  # Halogen rule
        X_flag = Xsingle and group_list[i] == 7 and period_list[i] <= 4

        ve_constr.add(group_list[i])
        prev_ve = group_list[i]  # previous valence electron

        chg_constr.add(t1[2 * i] - t1[2 * i + 1])
        lp_constr.add(t1[2 * i] - t1[2 * i + 1])
        ve_constr.add(-t1[2 * i] + t1[2 * i + 1])
        prev_ve += -t1_starts[2 * i] + t1_starts[2 * i + 1]

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
            # if X_constr is not None:
            #    X_constr.add(bo)

            prev_ve += (
                1
                + db_starts[bond_mapping[(a, b)]]
                + 2 * tb_starts[bond_mapping[(a, b)]]
            )

            if Xsingle and group_list[i] == 7 and period_list[i] <= 4:
                model.addConstr(bo == 1, name=f"XC_{i}")

            if i in M_list:
                model.addConstr(bo == 1, name=f"SB_{i}_{j}")

            # freeze the bond order for already octet atoms(period==2)
            # if bool(alreadyOctet[i]):
            #    # model.addConstr(eve_constr + group_list[i] == 8, name=f"eve_{i}")
            #    db[bond_mapping[(a, b)]].lb = db[bond_mapping[(a, b)]].Start
            #    db[bond_mapping[(a, b)]].ub = db[bond_mapping[(a, b)]].Start
            #    tb[bond_mapping[(a, b)]].lb = tb[bond_mapping[(a, b)]].Start
            #    tb[bond_mapping[(a, b)]].ub = tb[bond_mapping[(a, b)]].Start

        model.addConstr(lp_constr <= group_list[i], name=f"lp_{i}")

        # TODO: octet rule
        # if charged and period > 2, do not apply constraint
        # else, freeze the valence (octet rule)
        if not bool(overcharged[i]):
            # model.addConstr(ve_constr == prev_ve, name=f"ve_freeze_{i}") # don't know why this is not working
            model.addConstr(
                ve_constr == prev_ve, name=f"ve_freeze_{i}"
            )  # the same constraint with the maximize_bo
        else:
            model.addConstr(ve_constr >= prev_ve, name=f"ve_expanded_{i}")

        # freeze the formal charge also
        # if bool(alreadyOctet[i]):
        #    # model.addConstr(eve_constr + group_list[i] == 8, name=f"eve_{i}")
        #    t1[2 * i].lb = t1[2 * i].Start
        #    t1[2 * i].ub = t1[2 * i].Start
        #    t1[2 * i + 1].lb = t1[2 * i + 1].Start
        #    t1[2 * i + 1].ub = t1[2 * i + 1].Start

        # else:
        #    # apply loosen octet rule
        #    model.addConstr(eve_constr <= eve_list[i] - group_list[i], name=f"eve_{i}")
        #    # model.addConstr(eve_constr == eve_list[i] - group_list[i], name=f"eve_{i}")

        if eIsEven:
            # the number of valence electron is even number (no radical rule!)
            model.addConstr(ve_constr == 2 * even[i], name=f"noRad_{i}")

        # Ring Constraint
        if cleanUp and (i in ring_neighbors_info):
            for n, neighbor_bond_list in enumerate(ring_neighbors_info[i]):
                ring_constr = LinExpr()  # ring constraint
                for ring_bond in neighbor_bond_list:
                    ring_constr.add(
                        db[bond_mapping[ring_bond]] + tb[bond_mapping[ring_bond]]
                    )
                model.addConstr(ring_constr <= 1, name=f"ring_{i}_{n}")

        # halogen atoms have only one single bond
        # halogens might have any bond (halogen anions), and in such case, does not apply the constraint
        # if X_constr is not None and X_constr.size() > 0:
        #    model.addConstr(X_constr == 1, name=f"X_{i}")

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

    bo_priority = 2  # bond order maximization priority
    chg_priority = 3  # charge separation priority
    en_priority = 1  # electronegativity priority

    # if kwargs.get("mode", "") == "fc":
    #    bo_priority, chg_priority = chg_priority, bo_priority

    # model.setObjectiveN(
    #    max_bo_obj, 1, priority=bo_priority, weight=GRB.MAXIMIZE, name="max_bo"
    # )
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
        model.write(f"record_chg{stepIdx}.lp")

    model.optimize()

    # error handling
    if model.status != GRB.Status.OPTIMAL:
        print(
            f"resolve_chg: Optimization failed. (status: {model.status})",
            file=sys.stderr,
        )
        return np.zeros(atom_num), {}, (None, None, None)

    # result record
    if verbose:
        nSol = model.SolCount
        # print("nSol", nSol)
        for s in range(nSol):
            model.params.SolutionNumber = s
            model.write(f"output_chg{stepIdx}_{s}.sol")

    # retrieval
    bo_dict = {}
    chg_list = np.zeros(atom_num, dtype=np.int64)
    for i in range(bond_num):
        bo = 1 + int(db[i].X) + 2 * int(tb[i].X)
        bo_dict[bond_list[i]] = int(bo)
    for i in range(atom_num):
        chg_list[i] = int(t1[2 * i].X) - int(t1[2 * i + 1].X)

    db_values = [int(db[i].X) for i in range(bond_num)]
    tb_values = [int(tb[i].X) for i in range(bond_num)]
    t1_values = [int(t1[i].X) for i in range(2 * atom_num)]

    model.close()
    env.close()

    return chg_list, bo_dict, (db_values, tb_values, t1_values)


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
        ring_neighbors_info,
        en_list,
    ) = get_lists(molecule)

    atom_num, bond_num = len(z_list), len(bond_list)
    eIsEven = int(np.sum(z_list) - chg_mol) % 2 == 0
    resolve_step = 0
    kwargs["cleanUp"] = cleanUp

    chg_list0, bo_dict0, raw_outputs = maximize_bo(
        atom_num,
        bond_num,
        period_list,
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
    # print("Debug: chg_list0", chg_list0)
    # print("Debug: bo_dict0", bo_dict0)

    # early stop
    if chg_list0 is None and bo_dict0 is None:
        chg_list0, bo_matrix0 = None, None
    else:
        bo_matrix0 = np.zeros((atom_num, atom_num))
        for p, q in bo_dict0.keys():
            bo_matrix0[p][q] = bo_dict0[(p, q)]
            bo_matrix0[q][p] = bo_dict0[(p, q)]

    # check charge separation
    chg_sep = np.any(chg_list0 > 0) and np.any(chg_list0 < 0)

    bo_matrix1, chg_list1 = np.copy(bo_matrix0), np.copy(chg_list0)  # place holder

    # charge resolution
    if resolve and chg_sep:
        print("Debug: resolution")
        bo_sum = np.zeros(atom_num)
        for p, q in bo_dict0.keys():
            bo_sum[p] += bo_dict0[(p, q)]
            bo_sum[q] += bo_dict0[(p, q)]

        # alreadyOctet = (group_list + bo_sum - chg_list0 == 8) & (period_list == 2)
        overcharged = (period_list > 2) & (np.abs(chg_list0) != 0)
        print("Debug: overcharged", np.nonzero(overcharged))

        # ve_list is now expanded Valence list
        # eve_list = get_expanded_list(period_list, ve_list, chg_list, ring_list)
        eve_list = get_expanded_ve_list(period_list, group_list, ve_list, chg_list0)
        chg_list1, bo_dict1, raw_outputs1 = resolve_chg(
            atom_num,
            bond_num,
            period_list,
            group_list,
            bond_list,
            bond_mapping,
            eve_list,
            neighbor_list,
            en_list,
            ring_neighbors_info,
            chg_mol,
            eIsEven,
            overcharged,
            raw_outputs[0],
            raw_outputs[1],
            raw_outputs[2],
            stepIdx=resolve_step,
            **kwargs,
        )

        resolve_step += 1

        # error handling
        # if len(bo_dict) == 0:
        #    return None, None

        # if cleanUp:
        #    mve_list = get_modified_list(
        #        period_list, ve_list, chg_list, bo_dict, ring_list, ring_bond_list
        #    )
        #    bo_sum = np.zeros(atom_num)
        #    for p, q in bo_dict.keys():
        #        bo_sum[p] += bo_dict[(p, q)]
        #        bo_sum[q] += bo_dict[(p, q)]
        #    alreadyOctet = (group_list + bo_sum - chg_list == 8) & (period_list == 2)
        #    chg_list, bo_dict = resolve_chg(
        #        atom_num,
        #        bond_num,
        #        group_list,
        #        bond_list,
        #        bond_mapping,
        #        mve_list,
        #        neighbor_list,
        #        chg_mol,
        #        eIsEven,
        #        alreadyOctet,
        #        resolve_step,
        #        **kwargs,
        #    )
        #    resolve_step += 1

        # error handling
        if bo_dict1 is None and chg_list1 is None:
            chg_list1, bo_matrix1 = None, None
        else:
            bo_matrix1 = np.zeros((atom_num, atom_num))
            for p, q in bo_dict1.keys():
                bo_matrix1[p][q] = bo_dict1[(p, q)]
                bo_matrix1[q][p] = bo_dict1[(p, q)]

    return chg_list0, bo_matrix0, chg_list1, bo_matrix1


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

    chg_list, bo_dict, raw_outputs = maximize_bo(
        atom_num,
        bond_num,
        period_list,
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
    if bo_dict is None and chg_list is None:
        return None, None

    # check charge separation
    chg_sep = np.any(chg_list > 0) and np.any(chg_list < 0)

    # breakpoint()
    # charge resolution
    if resolve and chg_sep:
        # breakpoint()
        bo_sum = np.zeros(atom_num)
        for p, q in bo_dict.keys():
            bo_sum[p] += bo_dict[(p, q)]
            bo_sum[q] += bo_dict[(p, q)]

        # TODO: Check the condition for overcharged atoms
        # 1. period > 2
        # 2. non-zero charge on itself

        overcharged = (period_list > 2) & (np.abs(chg_list) != 0)
        # print("overcharged", np.nonzero(overcharged))
        # print("alreadyOctet", np.nonzero(alreadyOctet))

        # alreadyOctet = \
        #    (group_list + bo_sum - chg_list == 8) & \
        #    (period_list == 2) & \
        #    (np.abs(chg_list) == 0)

        # ve_list is now expanded Valence list
        # eve_list = get_expanded_list(period_list, ve_list, chg_list, ring_list)
        eve_list = get_expanded_ve_list(period_list, group_list, ve_list, chg_list)
        chg_list, bo_dict, raw_outputs2 = resolve_chg(
            atom_num,
            bond_num,
            period_list,
            group_list,
            bond_list,
            bond_mapping,
            eve_list,
            neighbor_list,
            en_list,
            ring_neighbors_info,
            chg_mol,
            eIsEven,
            overcharged,
            raw_outputs[0],
            raw_outputs[1],
            raw_outputs[2],
            stepIdx=resolve_step,
            **kwargs,
        )

        resolve_step += 1

        # error handling
        # if len(bo_dict) == 0:
        #    return None, None

        # if cleanUp:
        #    mve_list = get_modified_list(
        #        period_list, ve_list, chg_list, bo_dict, ring_list, ring_bond_list
        #    )
        #    bo_sum = np.zeros(atom_num)
        #    for p, q in bo_dict.keys():
        #        bo_sum[p] += bo_dict[(p, q)]
        #        bo_sum[q] += bo_dict[(p, q)]
        #    alreadyOctet = (group_list + bo_sum - chg_list == 8) & (period_list == 2)
        #    chg_list, bo_dict = resolve_chg(
        #        atom_num,
        #        bond_num,
        #        group_list,
        #        bond_list,
        #        bond_mapping,
        #        mve_list,
        #        neighbor_list,
        #        chg_mol,
        #        eIsEven,
        #        alreadyOctet,
        #        resolve_step,
        #        **kwargs,
        #    )
        #    resolve_step += 1

        # error handling
        if bo_dict is None and chg_list is None:
            return None, None

    bo_matrix = np.zeros((atom_num, atom_num))
    for p, q in bo_dict.keys():
        bo_matrix[p][q] = bo_dict[(p, q)]
        bo_matrix[q][p] = bo_dict[(p, q)]

    bo_listlist = bo_matrix.tolist()
    # breakpoint()

    return chg_list, bo_matrix


if __name__ == "__main__":
    import sys

    smi = sys.argv[1]
