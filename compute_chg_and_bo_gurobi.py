from gurobipy import *
from acerxn import chem
import numpy as np


def get_lists(molecule):
    period_list, group_list = molecule.get_period_group_list()
    adj_matrix = np.copy(molecule.get_matrix('adj'))
    adj_list = np.sum(adj_matrix, axis=1)
    ve_list = np.zeros_like(group_list)
    z_list = molecule.get_z_list()
    for i in range(len(group_list)):
        if period_list[i] == 1:
            ve_list[i] = group_list[i]
        else:
            ve_list[i] = 4
        if adj_list[i] > ve_list[i]:
            ve_list[i] = group_list[i]
    return period_list, group_list, z_list, ve_list, adj_list, molecule.get_bond_list(False)


def maximize_bo(molecule, chg_mol, mode):
    period_list, group_list, z_list, ve_list, adj_list, bond_list = get_lists(molecule)
    atom_num, bond_num = len(z_list), len(bond_list)

    if atom_num == 1:
        return np.array([chg_mol]), {}
    multiplicity = (np.sum(z_list)-chg_mol)%2+1

    model = Model('maximize_bo')
    bo = model.addVars(bond_num, name='bo', vtype=GRB.INTEGER)
    #chg = model.addVars(atom_num, lb=-np.inf, name='chg', vtype=GRB.INTEGER)
    t1 = model.addVars(2*atom_num, name='t1', vtype=GRB.INTEGER)
    if mode == 'heuristics':
        t2 = model.addVars(2*atom_num, name='t2', vtype=GRB.CONTINUOUS)
    if multiplicity == 1:
        lp = model.addVars(atom_num, name='lp', vtype=GRB.INTEGER)
 
    for i in range(bond_num):
        p, q = bond_list[i]
        model.addConstr(bo[i] <= min([ve_list[p]-adj_list[p], ve_list[q]-adj_list[q], 2]))

    constr = LinExpr()
    for i in range(atom_num):
        constr1 = LinExpr()
        constr2 = LinExpr(t1[2*i+1]-t1[2*i])
        constr3 = LinExpr(t1[2*i]-t1[2*i+1])
        if multiplicity == 1:
            constr4 = LinExpr(t1[2*i]-t1[2*i+1]+2*lp[i])
        constr.add(t1[2*i]-t1[2*i+1])
        for j in range(bond_num):
            if i in bond_list[j]:
                constr1.add(bo[j])
                constr2.add(bo[j])
                constr3.add(bo[j])    
                if multiplicity == 1:
                    constr4.add(bo[j])
        #model.addConstr(t1[2*i]-t1[2*i+1] == chg[i])
        if mode == 'heuristics':
            model.addConstr(t2[2*i]-t2[2*i+1] == (t1[2*i]-t1[2*i+1])+0.1*group_list[i]-0.4)
        model.addConstr(constr2 <= 2*ve_list[i]-adj_list[i]-group_list[i])
        model.addConstr(constr3 <= group_list[i]-adj_list[i])
        if multiplicity == 1:
            model.addConstr(constr4 == group_list[i]-adj_list[i])
    model.addConstr(constr == chg_mol)

    objective1 = LinExpr()
    objective2 = LinExpr()
    objective3 = LinExpr()
    for i in range(bond_num):
        objective1.add(bo[i])
    for i in range(atom_num):
        objective2.add(t1[2*i]+t1[2*i+1])
        if mode == 'heuristics':
            objective3.add((0.1*group_list[i]+0.6)*(t2[2*i]+t2[2*i+1]))
    model.setObjectiveN(objective1, 0, 2, GRB.MAXIMIZE)
    model.setObjectiveN(objective2, 1, 1, GRB.MINIMIZE)
    if mode == 'heuristics':
        model.setObjectiveN(objective3, 2, 0, GRB.MINIMIZE)
    model.setParam(GRB.Param.OutputFlag, 0)
    model.setParam(GRB.Param.TimeLimit, 1)
    model.optimize()

    if model.status != GRB.Status.OPTIMAL:
        return None, None

    bo_dict = {}
    chg_list = np.zeros(atom_num)
    for i in range(bond_num):
        bo_dict[bond_list[i]] = int(model.getVarByName('bo[{0}]'.format(i)).Xn)
    for i in range(atom_num):
        chg_list[i] = int(model.getVarByName('t1[{0}]'.format(2*i)).Xn)-int(model.getVarByName('t1[{0}]'.format(2*i+1)).Xn)
    
    return chg_list, bo_dict
 

def resolve_chg(molecule, chg_list, bo_dict, chg_mol, mode):
    period_list, group_list, z_list, ve_list, adj_list, bond_list = get_lists(molecule)
    atom_num, bond_num = len(z_list), len(bond_list)

    if atom_num == 1:
        return np.array([chg_mol]), {}

    multiplicity = (np.sum(z_list)-chg_mol)%2+1
    
    for i in range(atom_num):
        if (period_list[i] > 2 and group_list[i] > 4 and chg_list[i] > 0) or (group_list[i] < 4 and chg_list[i] < 0):
            ve_list[i] += chg_list[i]

    model = Model('resolve_chg')
    bo = model.addVars(bond_num, name='bo', vtype=GRB.INTEGER)
    #chg = model.addVars(atom_num, lb=-np.inf, name='chg', vtype=GRB.INTEGER)
    t1 = model.addVars(2*atom_num, name='t1', vtype=GRB.INTEGER)
    if mode == 'heuristics':
        t2 = model.addVars(2*atom_num, name='t2', vtype=GRB.CONTINUOUS)
    octet = model.addVars(atom_num, name='octet', vtype=GRB.INTEGER)
    if multiplicity == 1:
        lp = model.addVars(atom_num, name='lp', vtype=GRB.INTEGER)

    for i in range(bond_num):
        p, q = bond_list[i]
        model.addConstr(bo[i] <= min([ve_list[p]-adj_list[p], ve_list[q]-adj_list[q], 2]))

    constr = LinExpr()
    for i in range(atom_num):
        constr1 = LinExpr()
        constr2 = LinExpr(t1[2*i+1]-t1[2*i])
        constr3 = LinExpr(t1[2*i]-t1[2*i+1])
        if multiplicity == 1:
            constr4 = LinExpr(t1[2*i]-t1[2*i+1]+2*lp[i])
        #constr5 = LinExpr(t1[2*i+1]-t1[2*i]+octet[i])
        constr.add(t1[2*i]-t1[2*i+1])
        for j in range(bond_num):
            if i in bond_list[j]:
                constr1.add(bo[j])
                constr2.add(bo[j])
                constr3.add(bo[j])
                if multiplicity == 1:
                    constr4.add(bo[j])
                #constr5.add(bo[j])
        #model.addConstr(t1[2*i]-t1[2*i+1] == chg[i])
        if mode == 'heuristics':
            model.addConstr(t2[2*i]-t2[2*i+1] == (t1[2*i]-t1[2*i+1])+0.1*group_list[i]-0.4)
        model.addConstr(constr2 <= (2*ve_list[i]-adj_list[i]-group_list[i]))
        model.addConstr(constr3 <= group_list[i]-adj_list[i])
        if multiplicity == 1:
            model.addConstr(constr4 == group_list[i]-adj_list[i])
        #if ve_list[i] <= 4:
        #    model.addConstr(constr5 == (2*ve_list[i]-adj_list[i]-group_list[i]))
    model.addConstr(constr == chg_mol)

    objective1 = LinExpr()
    objective2 = LinExpr()
    objective3 = LinExpr()
    for i in range(atom_num):
        #objective1.add((2*ve_list[i]-group_list[i]-(t1[2*i]-t1[2*i+1])-3*adj_list[i])/2)
        objective2.add(3*(t1[2*i]+t1[2*i+1]))
        if mode == 'heuristics':
            objective3.add((0.1*group_list[i]+0.6)*(t2[2*i]+t2[2*i+1]))
    for i in range(bond_num):
        #objective1.add(-3*bo[i])
        objective2.add(-bo[i])
    #model.setObjectiveN(objective1, 0, 2, GRB.MINIMIZE)
    model.setObjectiveN(objective2, 1, 1, GRB.MINIMIZE)
    if mode == 'heuristics':
        model.setObjectiveN(objective3, 2, 0, GRB.MINIMIZE)
    model.setParam(GRB.Param.OutputFlag, 0)
    model.setParam(GRB.Param.TimeLimit, 1)
    #model.write('record.lp')
    model.optimize()

    if model.status != GRB.Status.OPTIMAL:
        return None, None

    bo_dict = {}
    chg_list = np.zeros(atom_num)
    lp_list = np.zeros(atom_num)
    for i in range(bond_num):
        bo_dict[bond_list[i]] = int(model.getVarByName('bo[{0}]'.format(i)).Xn)
    for i in range(atom_num):
        chg_list[i] = int(model.getVarByName('t1[{0}]'.format(2*i)).Xn)-int(model.getVarByName('t1[{0}]'.format(2*i+1)).Xn)
    for i in range(atom_num):
        lp_list[i] = int(model.getVarByName('lp[{0}]'.format(i)).Xn)

    return chg_list, bo_dict


def compute_chg_and_bo(molecule, chg_mol, resolve=True, mode='heuristics'):
    bo_matrix = np.copy(molecule.get_adj_matrix())
    chg_list, bo_dict = maximize_bo(molecule, chg_mol, mode)
    if not chg_list is None:
        if resolve:
            chg_list, bo_dict = resolve_chg(molecule, chg_list, bo_dict, chg_mol, mode)
        if not chg_list is None:
            for p, q in bo_dict.keys():
                bo_matrix[p][q] += bo_dict[(p, q)]
                bo_matrix[q][p] += bo_dict[(p, q)]
    
    return chg_list, bo_matrix
