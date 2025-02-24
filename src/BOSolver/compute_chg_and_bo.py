import sys

import pulp as pl
import numpy as np

from BOSolver.utils.chem import Chem, get_lists


def moSolve(prob, objs, verbose: bool, log_prefix=None):
    """multi-objective optimization. Simple modification of pulp.LpProblem.sequentialSolve"""
    statuses = []
    objvalues = []
    if not (prob.solver):
        prob.solver = pl.LpSolverDefault
        if not verbose:
            prob.solver.msg = False  # suppress the output
    for i, (_, obj, s) in enumerate(objs):
        prob.setObjective(obj)
        prob.sense = s
        status = prob.solver.actualSolve(prob)
        statuses.append(status)
        objvalues.append(obj.value())
        if verbose:
            prob.writeLP(f"{log_prefix}_record_obj{i}.lp")
        if s == pl.const.LpMinimize:
            prob += obj <= obj.value(), f"Obj_{i}"
        elif s == pl.const.LpMaximize:
            prob += obj >= obj.value(), f"Obj_{i}"
    return prob, statuses, objvalues


# FIXME: BOSolver cannot solve a system with odd number of electrons
def optimize_bo(
    atom_num,
    bond_num,
    period_list,
    ve_list,
    bond_list,
    bond_mapping,
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
    prob = pl.LpProblem("optimize_bo", pl.LpMaximize)

    ### option parsing
    verbose = kwargs.get("printOptLog", False)

    cleanUp = kwargs.get("cleanUp", False)
    RingConstr = cleanUp and (len(ring_neighbors_info) > 0)
    Xsingle = cleanUp and kwargs.get("HalogenConstraint", False)
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

    db = pl.LpVariable.dicts("dbFlag", range(bond_num), 0, 1, pl.LpBinary)
    tb = pl.LpVariable.dicts("tbFlag", range(bond_num), 0, 1, pl.LpBinary)

    prob.extend({f"BondOrderFlag_{i}": db[i] + tb[i] <= 1 for i in range(bond_num)})

    # t1: formal charge
    # t1[2i]: fc+ | t1[2i+1]: fc-
    # t1[2i] - t1[2i+1] : fc of atom i
    # t1[2i] + t1[2i+1] : abs(fc) of atom i
    t1 = pl.LpVariable.dicts("t1", range(2 * atom_num), 0, None, pl.LpInteger)

    # t2: formal charge for weighted objective function
    # weight considering electronegativity

    # o: octet distance
    # the distance between the number of valence elctrons and
    # the octet number(2 for 1st period, 8 for 2nd or higher period)
    # o = 8 - (g - c + b)

    # Set Initial Values
    for i in range(bond_num):
        db[i].setInitialValue(db_starts[i])
        tb[i].setInitialValue(tb_starts[i])
    for i in range(2 * atom_num):
        t1[i].setInitialValue(t1_starts[i])

    # TODO: Revision of Halogen Constraint
    # Halogen atoms, especially Cl and Br, are not allowed for
    # following the extended octet rule.
    # RDKit does not allow Cl and Br to have valence state greater than 1

    # even: dummy variable to force no. of electrons even
    even = pl.LpVariable.dicts("even", range(atom_num), 0, None, pl.LpInteger)

    ### objectives and constraints construction
    # objective functions
    min_od_obj = pl.LpAffineExpression(name="min_od")
    min_fc_obj = pl.LpAffineExpression(name="min_fc")
    max_bo_obj = pl.LpAffineExpression(name="max_bo")
    min_en_obj = pl.LpAffineExpression(name="min_en")

    # constraints
    chg_constr = pl.LpAffineExpression(name="chg_consv")

    for i in range(atom_num):
        lp_constr = pl.LpAffineExpression(name=f"lp_{i}")
        ve_constr = pl.LpAffineExpression(name=f"ve_{i}")

        ve_constr.addInPlace(ve_list[i])

        chg_constr.addInPlace(t1[2 * i] - t1[2 * i + 1])
        lp_constr.addInPlace(t1[2 * i] - t1[2 * i + 1])
        ve_constr.addInPlace(-t1[2 * i] + t1[2 * i + 1])
        min_fc_obj.addInPlace(t1[2 * i] + t1[2 * i + 1])
        min_en_obj.addInPlace(en_list[i] * (t1[2 * i] - t1[2 * i + 1]))

        # summation over bond
        for j in neighbor_list[i]:
            a, b = i, j
            if a > b:
                a, b = b, a

            bo = pl.LpAffineExpression(
                1 + db[bond_mapping[(a, b)]] + tb[bond_mapping[(a, b)]] * 2
            )

            lp_constr.addInPlace(bo)
            ve_constr.addInPlace(bo)

            max_bo_obj.addInPlace(bo)

            # halogen atoms have only one single bond
            # halogens might have any bond (halogen anions), and in such case, does not apply the constraint
            if Xsingle and ve_list[i] == 7 and period_list[i] <= 4:
                prob += bo == 1, f"XC_{i}"

            # metal constraint
            if i in M_list:
                prob += bo == 1, f"SB_{i}_{j}"

        # the number of lone pair should not be negative
        prob += lp_constr <= ve_list[i], f"lp_{i}"

        # octet rule
        # octet distance
        if period_list[i] == 1:
            min_od_obj.addInPlace(2 - ve_constr)
            prob += 2 - ve_constr >= 0, f"od_{i}"
        elif period_list[i] == 2 or len(neighbor_list[i]) <= 4:
            min_od_obj.addInPlace(8 - ve_constr)
            prob += 8 - ve_constr >= 0, f"od_{i}"

        # the number of valence electron is even number (no radical rule!)
        if eIsEven:
            prob += ve_constr == 2 * even[i], f"noRad_{i}"

        # Ring Constraint
        if RingConstr and (i in ring_neighbors_info):
            for n, neighbor_bond_list in enumerate(ring_neighbors_info[i]):
                ring_constr = pl.LpAffineExpression(name=f"ring_{i}_{n}")
                for ring_bond in neighbor_bond_list:
                    ring_constr.addInPlace(
                        db[bond_mapping[ring_bond]] + tb[bond_mapping[ring_bond]]
                    )
                prob += ring_constr <= 1, f"ring_{i}_{n}"

    prob += chg_constr == chg_mol, "chg_consv"

    od_priority = 1  # octet distance priority
    chg_priority = 2  # charge separation priority
    bo_priority = 3  # bond order maximization priority
    en_priority = 4  # electronegativity priority

    if kwargs.get("fcmode", False):
        od_priority, chg_priority = chg_priority, od_priority

    # TODO: use electron negativity as the objective function on the request
    objs = [
        (od_priority, min_od_obj, pl.LpMinimize),
        (chg_priority, min_fc_obj, pl.LpMinimize),
        (bo_priority, max_bo_obj, pl.LpMaximize),
        # (en_priority, min_en_obj, pl.LpMinimize),
    ]
    objs = sorted(objs, key=lambda x: x[0])

    # Pulp optimization
    prob, statuses, objvalues = moSolve(prob, objs, verbose, "optimize_bo")

    # error handling
    for i, status in enumerate(statuses):
        if status != pl.LpStatusOptimal:
            print(
                f"optimize_bo: Obj{i} Optimization failed. (status: {pl.LpStatus[status]})",
                file=sys.stderr,
            )
            return None, None, (None, None, None)

    # result record
    if verbose:
        import json

        output = prob.toDict()
        output["status"] = statuses
        output["obj_values"] = objvalues
        json.dump(output, open("output.json", "w"), indent=4, default=str)

    # retrieval
    bo_dict = {}
    chg_list = np.zeros(atom_num, dtype=np.int64)
    for i in range(bond_num):
        bo = 1 + db[i].value() + 2 * tb[i].value()
        bo_dict[bond_list[i]] = int(bo)
    for i in range(atom_num):
        chg_list[i] = int(t1[2 * i].value() - t1[2 * i + 1].value())
    db_values = [int(db[i].value()) for i in range(bond_num)]
    tb_values = [int(tb[i].value()) for i in range(bond_num)]
    t1_values = [int(t1[i].value()) for i in range(2 * atom_num)]

    return chg_list, bo_dict, (db_values, tb_values, t1_values)


def resolve_chg(
    atom_num,
    bond_num,
    period_list,
    ve_list,
    bond_list,
    bond_mapping,
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
    prob = pl.LpProblem(f"resolve_chg{stepIdx}", pl.LpMaximize)

    verbose = kwargs.get("printOptLog", False)
    Xsingle = kwargs.get("HalogenConstraint", False)
    cleanUp = kwargs.get("cleanUp", False)
    RingConstr = cleanUp and (len(ring_neighbors_info) > 0)
    M_list = kwargs.get("MetalCenters", [])

    # bo: bond order
    db = pl.LpVariable.dicts("dbFlag", range(bond_num), 0, 1, pl.LpBinary)
    tb = pl.LpVariable.dicts("tbFlag", range(bond_num), 0, 1, pl.LpBinary)
    prob.extend({f"BondOrderFlag_{i}": db[i] + tb[i] <= 1 for i in range(bond_num)})

    # t1: formal charge
    t1 = pl.LpVariable.dicts("t1", range(2 * atom_num), 0, None, pl.LpInteger)

    # t2: formal charge for weighted objective function
    # weight considering electronegativity
    # t2 = model.addVars(2 * atom_num, name="t2", vtype=GRB.CONTINUOUS)

    # Set Initial Values
    for i in range(bond_num):
        db[i].setInitialValue(db_starts[i])
        tb[i].setInitialValue(tb_starts[i])
    for i in range(2 * atom_num):
        t1[i].setInitialValue(t1_starts[i])

    # even: dummy variable to force no. of electrons even
    even = pl.LpVariable.dicts("even", range(atom_num), 0, None, pl.LpInteger)

    ### objectives and constraints construction
    # objective functions
    min_fc_obj = pl.LpAffineExpression(name="min_fc")
    max_bo_obj = pl.LpAffineExpression(name="max_bo")
    min_en_obj = pl.LpAffineExpression(name="min_en")
    # constraints
    chg_constr = pl.LpAffineExpression(name="chg_consv")

    for i in range(atom_num):
        lp_constr = pl.LpAffineExpression(name=f"lp_{i}")
        ve_constr = pl.LpAffineExpression(name=f"ve_{i}")
        X_flag = Xsingle and ve_list[i] == 7 and period_list[i] <= 4

        ve_constr.addInPlace(ve_list[i])
        prev_ve = ve_list[i]  # previous valence electron

        chg_constr.addInPlace(t1[2 * i] - t1[2 * i + 1])
        lp_constr.addInPlace(t1[2 * i] - t1[2 * i + 1])
        ve_constr.addInPlace(-t1[2 * i] + t1[2 * i + 1])
        min_fc_obj.addInPlace(t1[2 * i] + t1[2 * i + 1])
        min_en_obj.addInPlace(en_list[i] * (t1[2 * i] - t1[2 * i + 1]))
        prev_ve += -t1_starts[2 * i] + t1_starts[2 * i + 1]  # previous valence electron

        # summation over bond
        for j in neighbor_list[i]:
            a, b = i, j
            if a > b:
                a, b = b, a

            bo = pl.LpAffineExpression(
                1 + db[bond_mapping[(a, b)]] + tb[bond_mapping[(a, b)]] * 2
            )

            lp_constr.addInPlace(bo)
            ve_constr.addInPlace(bo)

            max_bo_obj.addInPlace(bo)

            prev_ve += (
                1
                + db_starts[bond_mapping[(a, b)]]
                + 2 * tb_starts[bond_mapping[(a, b)]]
            )  # previous valence electron
            # Halogen Constraint
            # halogen atoms should obey the octet rule
            # (no extended octet rule for halogens)
            # TODO: Revision of Halogen Constraint
            # Halogen atoms, especially Cl and Br, are not allowed for
            # following the extended octet rule.
            # RDKit does not allow Cl and Br to have valence state greater than 1

            # halogen atoms have only one single bond
            # halogens might not have any bond (halogen anions), and in such case, does not apply the constraint
            if Xsingle and ve_list[i] == 7 and period_list[i] <= 4:
                prob += bo == 1, f"XC_{i}"

            # metal constraint
            if i in M_list:
                prob += bo == 1, f"SB_{i}_{j}"

        # the number of lone pair should not be negative
        prob += lp_constr <= ve_list[i], f"lp_{i}"

        # octet rule
        # if charged and period > 2, apply expanded octet rule
        # else, freeze the valence (octet rule)
        if not bool(overcharged[i]):
            prob += (
                ve_constr == prev_ve,
                f"ve_freeze_{i}",
            )  # don't know why this is not working
            # prob += ve_constr <= prev_ve, f"ve_freeze_{i}"  # this constraint goes wrong with azide moiety
        else:
            prob += ve_constr >= prev_ve, f"ve_expanded_{i}"

        # the number of valence electron is even number (no radical rule!)
        if eIsEven:
            prob += ve_constr == 2 * even[i], f"noRad_{i}"

        # Ring Constraint
        if RingConstr and (i in ring_neighbors_info):
            for n, neighbor_bond_list in enumerate(ring_neighbors_info[i]):
                ring_constr = pl.LpAffineExpression(name=f"ring_{i}_{n}")
                for ring_bond in neighbor_bond_list:
                    ring_constr.addInPlace(
                        db[bond_mapping[ring_bond]] + tb[bond_mapping[ring_bond]]
                    )
                prob += ring_constr <= 1, f"ring_{i}_{n}"

    prob += chg_constr == chg_mol, "chg_consv"

    ### optimization
    chg_priority = 1  # charge separation priority
    bo_priority = 2  # bond order maximization priority
    en_priority = 3  # electronegativity priority

    objs = [
        # (bo_priority, max_bo_obj, pl.LpMaximize),
        (chg_priority, min_fc_obj, pl.LpMinimize),
        (en_priority, min_en_obj, pl.LpMinimize),
    ]
    objs = sorted(objs, key=lambda x: x[0])

    # Pulp optimization
    prob, statuses, objvalues = moSolve(prob, objs, verbose, "resolve_chg")

    # error handling
    for i, status in enumerate(statuses):
        if status != pl.LpStatusOptimal:
            print(
                f"resolve_chg: Obj{i} Optimization failed. (status: {pl.LpStatus[status]})",
                file=sys.stderr,
            )
            return None, None, (None, None, None)

    # result record
    if verbose:
        import json

        output = prob.toDict()
        output["status"] = statuses
        output["obj_values"] = objvalues
        json.dump(
            output, open(f"output_resolve{stepIdx}.json", "w"), indent=4, default=str
        )

    # retrieval
    bo_dict = {}
    chg_list = np.zeros(atom_num, dtype=np.int64)
    for i in range(bond_num):
        bo = 1 + db[i].value() + 2 * tb[i].value()
        bo_dict[bond_list[i]] = int(bo)
    for i in range(atom_num):
        chg_list[i] = int(t1[2 * i].value() - t1[2 * i + 1].value())
    db_values = [int(db[i].value()) for i in range(bond_num)]
    tb_values = [int(tb[i].value()) for i in range(bond_num)]
    t1_values = [int(t1[i].value()) for i in range(2 * atom_num)]

    return chg_list, bo_dict, (db_values, tb_values, t1_values)


def compute_chg_and_bo(
    molecule: Chem.Mol, chg_mol, resolve=True, cleanUp=True, **kwargs
):
    """
    Compute the charge and bond order for a given molecule.

    Args:
        molecule (Chem.Mol): The RDKit Mol object containing atomic and bonding information (except bond orders).
        chg_mol (int): The total charge of the molecule.
        resolve (bool, optional): Whether to go through charge resolution step if needed. Defaults to True.
        cleanUp (bool, optional): Whether to apply heuristics that cleans up the resulting molecular graph. Defaults to True.
        **kwargs: Additional keyword arguments to be passed to the optimize_bo and resolve_chg functions.

        kwargs include:
        HaloGenConstraint (bool, optional): Whether to apply the halogen constraint (Halogens are considered to be terminal). Defaults to False.
                                            If cleanUp is False, this constraint is not applied.
        MetalCenters (list, optional): The list of atom indices that are considered to be metal centers.
                                    Bonds of metal atoms are forced to be single bonds. Defaults to [].
        fcmode (bool, optional): Whether to set formal charge separation minimization as the primary objective. Defaults to False.


    Returns:
        chg_list: the list of formal charges for each atom
        bo_matrix: bond order matrix
    """

    (
        period_list,
        ve_list,
        z_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        ring_neighbors_info,
        en_list,
    ) = get_lists(molecule)

    atom_num, bond_num = len(z_list), len(bond_list)
    eIsEven = int(np.sum(z_list) - chg_mol) % 2 == 0
    resolve_step = 0
    kwargs["cleanUp"] = cleanUp  # suppress

    chg_list, bo_dict, raw_outputs = optimize_bo(
        atom_num,
        bond_num,
        period_list,
        ve_list,
        bond_list,
        bond_mapping,
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

    # charge resolution
    if resolve and chg_sep:
        bo_sum = np.zeros(atom_num)
        for p, q in bo_dict.keys():
            bo_sum[p] += bo_dict[(p, q)]
            bo_sum[q] += bo_dict[(p, q)]

        # TODO: Check the condition for overcharged atoms
        # 1. period > 2
        # 2. non-zero charge on itself
        overcharged = (period_list > 2) & (np.abs(chg_list) != 0)

        chg_list, bo_dict, _ = resolve_chg(
            atom_num,
            bond_num,
            period_list,
            ve_list,
            bond_list,
            bond_mapping,
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

        ### TODO: Do resolve_chg several times if needed

        # error handling
        if bo_dict is None and chg_list is None:
            return None, None

    bo_matrix = np.zeros((atom_num, atom_num))

    assert bo_dict is not None
    assert chg_list is not None

    for p, q in bo_dict.keys():
        bo_matrix[p][q] = bo_dict[(p, q)]
        bo_matrix[q][p] = bo_dict[(p, q)]

    return chg_list, bo_matrix


def compute_chg_and_bo_debug(molecule, chg_mol, resolve=True, cleanUp=True, **kwargs):
    (
        period_list,
        ve_list,
        z_list,
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

    chg_list0, bo_dict0, raw_outputs = optimize_bo(
        atom_num,
        bond_num,
        period_list,
        ve_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        en_list,
        ring_neighbors_info,
        chg_mol,
        eIsEven,
        **kwargs,
    )

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

        # TODO: Check the condition for overcharged atoms
        # 1. period > 2
        # 2. non-zero charge on itself
        overcharged = (period_list > 2) & (np.abs(chg_list0) != 0)
        print("Debug: overcharged", np.nonzero(overcharged))

        chg_list1, bo_dict1, raw_outputs1 = resolve_chg(
            atom_num,
            bond_num,
            period_list,
            ve_list,
            bond_list,
            bond_mapping,
            neighbor_list,
            en_list,
            ring_neighbors_info,
            chg_mol,
            eIsEven,
            overcharged,
            raw_outputs[0],
            raw_outputs[1],
            raw_outputs[2],
            **kwargs,
        )

        resolve_step += 1

        # error handling
        if bo_dict1 is None and chg_list1 is None:
            chg_list1, bo_matrix1 = None, None
        else:
            bo_matrix1 = np.zeros((atom_num, atom_num))
            for p, q in bo_dict1.keys():
                bo_matrix1[p][q] = bo_dict1[(p, q)]
                bo_matrix1[q][p] = bo_dict1[(p, q)]

    return chg_list0, bo_matrix0, chg_list1, bo_matrix1
