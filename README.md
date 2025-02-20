```

                                                                                                                                                
BBBBBBBBBBBBBBBBB        OOOOOOOOO        SSSSSSSSSSSSSSS                  lllllll                                                              
B::::::::::::::::B     OO:::::::::OO    SS:::::::::::::::S                 l:::::l                                                              
B::::::BBBBBB:::::B  OO:::::::::::::OO S:::::SSSSSS::::::S                 l:::::l                                                              
BB:::::B     B:::::BO:::::::OOO:::::::OS:::::S     SSSSSSS                 l:::::l                                                              
  B::::B     B:::::BO::::::O   O::::::OS:::::S               ooooooooooo    l::::lvvvvvvv           vvvvvvv eeeeeeeeeeee    rrrrr   rrrrrrrrr   
  B::::B     B:::::BO:::::O     O:::::OS:::::S             oo:::::::::::oo  l::::l v:::::v         v:::::vee::::::::::::ee  r::::rrr:::::::::r  
  B::::BBBBBB:::::B O:::::O     O:::::O S::::SSSS         o:::::::::::::::o l::::l  v:::::v       v:::::ve::::::eeeee:::::eer:::::::::::::::::r 
  B:::::::::::::BB  O:::::O     O:::::O  SS::::::SSSSS    o:::::ooooo:::::o l::::l   v:::::v     v:::::ve::::::e     e:::::err::::::rrrrr::::::r
  B::::BBBBBB:::::B O:::::O     O:::::O    SSS::::::::SS  o::::o     o::::o l::::l    v:::::v   v:::::v e:::::::eeeee::::::e r:::::r     r:::::r
  B::::B     B:::::BO:::::O     O:::::O       SSSSSS::::S o::::o     o::::o l::::l     v:::::v v:::::v  e:::::::::::::::::e  r:::::r     rrrrrrr
  B::::B     B:::::BO:::::O     O:::::O            S:::::So::::o     o::::o l::::l      v:::::v:::::v   e::::::eeeeeeeeeee   r:::::r            
  B::::B     B:::::BO::::::O   O::::::O            S:::::So::::o     o::::o l::::l       v:::::::::v    e:::::::e            r:::::r            
BB:::::BBBBBB::::::BO:::::::OOO:::::::OSSSSSSS     S:::::So:::::ooooo:::::ol::::::l       v:::::::v     e::::::::e           r:::::r            
B:::::::::::::::::B  OO:::::::::::::OO S::::::SSSSSS:::::So:::::::::::::::ol::::::l        v:::::v       e::::::::eeeeeeee   r:::::r            
B::::::::::::::::B     OO:::::::::OO   S:::::::::::::::SS  oo:::::::::::oo l::::::l         v:::v         ee:::::::::::::e   r:::::r            
BBBBBBBBBBBBBBBBB        OOOOOOOOO      SSSSSSSSSSSSSSS      ooooooooooo   llllllll          vvv            eeeeeeeeeeeeee   rrrrrrr            
                                                                                                                                                
                                                                                                                                                
```

# Bond Order Solver utilizing linear integer programming

BOSolver calculates the CORRECT bond orders of bonds between atoms from XYZ file of molecule(s).
Lewis diagram drawing procedure is translated into Integer Linear Programming (ILP),
and BOSolver solves the bond order assignment problem EXACTLY when elements of atoms, connectivity, and total charge is provided.

# Installation

## Conda

Installation via conda is recommended.
`conda install -c conda-forge bosolver`

## Pip

Installation via pip (PyPi) is also available, although dependencies, including `rdkit`, might not be installed properly.

## From source

Installation from source is also possible. Clone the repository and install with pip.

```git clone https://github.com/DOCH-2/BOSolver.git```

and then run

```pip install .```

# Usage

BOsolver requires an .xyz file (or text formatted in xyz) of a molecule (system) and the total charge of the molecule (system).
BOSolver can be used as a command line tool or as a Python package.

## as Command Line Tool

```bash
bosolve molecule.xyz 0
```

or pass the content of .xyz file directly

```bash
bosolve "$(cat molecule.xyz)" 0

```

For more details, run `bosolve -h`

## as Python package

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

# Reference and Citation

This work is based on the Doctoral thesis of the Dr. Kyunghoon Lee, one of the authors.
To cite this work, cite as the following:
