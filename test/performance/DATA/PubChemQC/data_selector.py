"""Randomly choose 1M from the whole dataset"""
import pathlib
import pandas as pd
import numpy as np
from rdkit import Chem

np.random.seed(42)

# neutral
fname = sorted(pathlib.Path(".").glob("neutral_tot.woRadical*"))[-1]
df = pd.read_csv(str(fname))
neutral_sel = df.sample(n=100000, random_state=42)
neutral_sel.to_csv(str(fname).replace("tot", "sel"),index=False, header=False)

# df = pd.read_csv("name_smi_neutral_limited", names=["cid"])
# neutral_smi_sel = df.sample(n=100000, random_state=42)
# neutral_smi_sel.to_csv("neutral_label", index=False, header=False)
