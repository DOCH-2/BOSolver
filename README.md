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

## Reference

This work is based on the Doctoral thesis of the Dr. Kyunghoon Lee, one of the authors.
