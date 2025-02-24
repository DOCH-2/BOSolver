"""Randomly choose 1M from the whole dataset"""

import pandas as pd
import numpy as np
from rdkit import Chem

np.random.seed(42)

# neutral
df = pd.read_csv("neutral_tot.woRadical250221")
neutral_sel = df.sample(n=100000, random_state=42)
neutral_sel.to_csv("neutral_sel.woRadical250221", index=False, header=False)

# df = pd.read_csv("name_smi_neutral_limited", names=["cid"])
# neutral_smi_sel = df.sample(n=100000, random_state=42)
# neutral_smi_sel.to_csv("neutral_label", index=False, header=False)
