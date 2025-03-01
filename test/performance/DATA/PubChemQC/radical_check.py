import pathlib
from rdkit import Chem

PT = Chem.GetPeriodicTable()
# with open("./neutral_tot.limited.withRadical") as f:
fname = sorted(pathlib.Path(".").glob("neutral_tot.withRadical*"))[-1]
with open(fname) as f:
    for xyz in f:
        fdir = xyz.strip()
        symbols = []
        with open(fdir) as file:
            for line_number, line in enumerate(file):
                if line_number == 0:
                    continue
                elif line_number == 1:
                    continue
                else:
                    symbol, _, _, _ = line.split()
                    symbols.append(symbol)
        z_list = [PT.GetAtomicNumber(symbol) for symbol in symbols]
        if sum(z_list) % 2 == 0:
            print(fdir)

# with open("./anion_tot.limited.withRadical") as f:
#    for xyz in f:
#        fdir, chg = xyz.strip().split()
#        symbols = []
#        with open(fdir) as file:
#            for line_number, line in enumerate(file):
#                if line_number == 0:
#                    continue
#                elif line_number == 1:
#                    continue
#                else:
#                    symbol, _, _, _ = line.split()
#                    symbols.append(symbol)
#        z_list = [PT.GetAtomicNumber(symbol) for symbol in symbols]
#        if (sum(z_list) - int(chg)) % 2 == 1:
#            print(fdir)

# with open("./cation_tot.limited.withRadical") as f:
#    for xyz in f:
#        fdir, chg = xyz.strip().split()
#        symbols = []
#        with open(fdir) as file:
#            for line_number, line in enumerate(file):
#                if line_number == 0:
#                    continue
#                elif line_number == 1:
#                    continue
#                else:
#                    symbol, _, _, _ = line.split()
#                    symbols.append(symbol)
#        z_list = [PT.GetAtomicNumber(symbol) for symbol in symbols]
#        if (sum(z_list) - int(chg)) % 2 == 0:
#            print(fdir, chg)
