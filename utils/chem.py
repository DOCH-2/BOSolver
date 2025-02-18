from collections import namedtuple
import numpy as np
from rdkit import Chem


PT = Chem.GetPeriodicTable()

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
    # "P": 3.053, # phosphorus, adjusted
    "S": 2.589,
    # "S": 3.089, # sulfur, adjusted
    "Cl": 2.869,
    # "Cl": 3.369, # Halogen, adjusted
    "Ar": 3.242,
    "K": 0.734,
    "Ca": 1.034,
    "Ga": 1.756,
    "Ge": 1.994,
    "As": 2.211,
    "Se": 2.424,
    "Br": 2.685,
    # "Br": 3.285, # Halogen, adjusted
    "Kr": 2.966,
    "Rb": 0.706,
    "Sr": 0.963,
    "In": 1.656,
    "Sn": 1.824,
    "Sb": 1.984,
    "Te": 2.158,
    "I": 2.359 + 1,
    # "I": 3.159, # Halogen, adjusted
    "Xe": 2.582,
}  # J. Am. Chem. Soc. 1989, 111, 25, 9003–9014


class MolData(
    namedtuple(
        "MolData",
        [
            "period_list",
            "ve_list",
            "z_list",
            "bond_list",
            "bond_mapping",
            "neighbor_list",
            "ring_neighbors_info",
            "en_list",
        ],
    )
):
    def __new__(
        cls,
        period_list,
        ve_list,
        z_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        ring_neighbors_info,
        en_list,
    ):
        return super().__new__(
            cls,
            period_list,
            ve_list,
            z_list,
            bond_list,
            bond_mapping,
            neighbor_list,
            ring_neighbors_info,
            en_list,
        )


def get_period_ve_list(mol: Chem.Mol):
    period_list = np.array([PT.GetRow(atom.GetAtomicNum()) for atom in mol.GetAtoms()])
    ve_list = np.array(
        [PT.GetNOuterElecs(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
    )

    return period_list, ve_list


def get_bond_info(mol: Chem.Mol):
    neighbor_list = [
        list(map(lambda x: x.GetIdx(), atom.GetNeighbors())) for atom in mol.GetAtoms()
    ]
    rd_bonds = list(mol.GetBonds())
    bond_list = list(map(lambda b: (b.GetBeginAtomIdx(), b.GetEndAtomIdx()), rd_bonds))
    bond_mapping = {key: val for key, val in enumerate(bond_list)}

    return neighbor_list, bond_list, bond_mapping


def get_ring_info(mol: Chem.Mol):
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
    new_z_list = np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])
    degree = [atom.GetDegree() for atom in mol.GetAtoms()]
    new_z_list[degree == 1] = 1  # H
    new_z_list[degree == 2] = 8  # O
    new_z_list[degree == 3] = 7  # N
    new_z_list[degree == 4] = 6  # C
    new_z_list[degree == 5] = 15  # P
    new_z_list[degree == 6] = 16  # S

    new_rd = Chem.rdchem.RWMol()
    for z in new_z_list:
        new_rd.AddAtom(Chem.rdchem.Atom(int(z)))
    for b in mol.GetBonds():
        new_rd.AddBond(
            b.GetBeginAtomIdx(), b.GetEndAtomIdx(), Chem.rdchem.BondType.SINGLE
        )
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

    return atoms_in_ring, bonds_in_ring, ring_neighbors_info


def get_lists(mol: Chem.Mol):
    # period, group
    period_list, ve_list = get_period_ve_list(mol)

    # atomic number
    z_list = np.array([atom.GetAtomicNum() for atom in mol.GetAtoms()])

    # neighbor, bond, bond_mapping
    neighbor_list, bond_list, bond_mapping = get_bond_info(mol)

    # ring info
    _, _, ring_neighbors_info = get_ring_info(mol)

    # electronegativity
    en_list = np.array([EN_TABLE[PT.GetElementSymbol(int(z))] for z in z_list])

    return MolData(
        period_list,
        ve_list,
        z_list,
        bond_list,
        bond_mapping,
        neighbor_list,
        ring_neighbors_info,
        en_list,
    )
