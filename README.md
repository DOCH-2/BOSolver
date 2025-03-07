# Bond Order Solver utilizing Integer Linear Programming

BOSolver calculates the CORRECT bond orders of bonds between atoms from XYZ file
of molecule(s).

Lewis diagram drawing procedure is translated into Integer Linear Programming
(ILP), and BOSolver solves the bond order assignment problem EXACTLY
when elements of atoms, connectivity, and total charge is provided.

Thanks to the ILP formulation, BOSolver has multiple strength points such as,

- BOSolver can handle complex molecules with
multiple resonance structures within *definite* time.
- BOSolver preserves the total charge of the molecule.
- BOSolver is free from loopholes of heuristics-based methods.
- BOSolver shows, empirically, stable time complexity
regardless of cases (worst-case, best-case, and average).

BOSolver relies on `RDKit` for the molecular representation,
and BOSolver can assign correct bond orders to `RDKit.Chem.Mol` objects
within a few lines.
Read the usage section for more details.

## Installation

### Pip

Installation via pip (PyPi) is available.

```bash
>>> pip install BOSolver
```

### From source

Installation from source is also possible. Clone the repository and
install with pip.

```git clone https://github.com/DOCH-2/BOSolver.git```

and then run

```bash
>>> cd BOSolver # move to the top directory
>>> python -m build
>>> pip install dist/BOSolver.xxx.whl
```

To check whether the installation is properly done,
you can use `pytest`. Run `pytest` at the top directory of the repository.

```bash
>>> cd BOSolver # move to the top directory
>>> pytest
```

## Usage

BOsolver requires an .xyz file (or text formatted in xyz) of a molecule (system)
and the total charge of the molecule (system).

BOSolver can be used as a command line tool or as a Python package.

### as Command Line Tool

```bash
bosolve molecule.xyz 0
```

or pass the content of .xyz file directly

```bash
bosolve "$(cat molecule.xyz)" 0
```

For more details, run `bosolve -h`

### as Python package

To assign bond orders to a rdkit.Chem.Mol object, use `BOSolver.bosolve.assignBO`

```python
from BOsolver.bosolve import assignBO, perceiveConn

mol = Chem.MolFromXYZFile("molecule.xyz")
chg = 0

# if molecule has no connectivity information, call perceiveConn first
if not mol.GetNumBonds():
    mol = perceiveConn(mol)

# if molecule has no connectivity information, then assignBO will raise an error
mol = assignBO(mol, chg)
```

### Notes

- A note on `perceiveConn` function

Providing the correct connectivity information is crucial for BOSolver to work correctly.

`perceiveConn` is a function that perceives connectivity information from the atomic coordinates.

Various algorithms can be applied to perceive the connectivity information.

Here is an example of using `perceiveConn` function using different perception algorithms.

```python
from BOsolver.bosolve import perceiveConn
from BOSolver.utils.coord2adj import BondPerception, CovalentRadius, Simple, BaberHodgkin

mol = Chem.MolFromXYZFile("molecule.xyz")
chg = 0

algo1 = CovalentRadius(relTol=1.1, absTol=0.0) # default
mol = perceiveConn(mol, algo1) # apply CovalentRadius algorithm

algo2 = Simple() # relTol=1.15, absTol=0.0
mol = perceiveConn(mol, algo2) # apply Simple algorithm

algo3 = BaberHodgkin() # relTol=1.0, absTol=0.45
mol = perceiveConn(mol, algo3) # apply BaberHodgkin algorithm

# or define your own algorithm
# refer to the BondPerception class for more details
class MyAlgorithm(BondPerception):
    def __init__(self, *args, **kwargs):
        super().__init__()

    # write your algorithm
    # blah blah blah
    # blah blah blah

    def __call__(self, mol):
        # your own algorithm
        return mol
```

Class `CovalenRadius` and its subclasses has an attribute `R_TABLE`,
which contains covalent radii of elements.

As you are dealing with specific systems, you might want to adjust
the covalent radii of some elements.

```python
from BOSolver.utils.coord2adj import CovalentRadius

algo = CovalentRadius()
algo.R_TABLE[6] = 0.8 # change the covalent radius of Carbon(atomic number 6) to 0.8
```

- A note on dealing with a system with unknown charge

BOSolver requires the total charge of the system (it does not have to be single molecule)
to work correctly.

When the total charge of the system is unknown, you can try any charge values to guess the total
charge of the system.

In practice, for a given connectivity, there are only one or two possible charge values that
give a valid bond order assignment.
So try giving charge values from -2 to 2 to find the correct charge value.

```python
from rdkit import Chem
from BOsolver.bosolve import assignBO
chg_table = [0, 1, -1, 2, 2] # start from small absolute value
candidates = []

mol = Chem.MolFromXYZFile("molecule.xyz")
for chg in chg_table:
    try:
        mol = assignBO(mol, chg)
        Chem.SanitizeMol(mol) # RDkit's SanitizeMol function will check the validity
    except:
        continue
    else:
        candidates.append(mol)
```

## Reference

This work is based on the Doctoral thesis of the Dr. Kyunghoon Lee, one of the authors.

## Support

If you have any questions or suggestions, open issues or contact Jinwon Lee.
