"""analyze result of test.py"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import seaborn as sns


def append_ChemVal_ChemStab(df):
    df["ChemVal"] = df["NoRadical"] + df["Sanitize"] + df["ChgConserve"]
    df["ChemStab"] = df["BOSum"] - df["ChgSepa"] + df["nAromatic"]
    return df


def get_chemstab(df: pd.Series):
    new_df = pd.DataFrame(
        {
            "nOctet": df[
                ["nOctet_obb", "nOctet_igx", "nOctet_x2m", "nOctet_bos"]
            ].to_list(),
            "ChgSepa": df[
                ["ChgSepa_obb", "ChgSepa_igx", "ChgSepa_x2m", "ChgSepa_bos"]
            ].to_list(),
            "nAromatic": df[
                [
                    "nAromatic_obb",
                    "nAromatic_igx",
                    "nAromatic_x2m",
                    "nAromatic_bos",
                ]
            ].to_list(),
        },
        index=["OBB", "IGX", "X2M", "BOS"],
    )

    cols = np.array(["bos", "x2m", "igx", "obb"])

    # hierarchical comparision
    # 1. ChgSepa
    # 2. BOSum
    # 3. nAromatic
    # 4. nOctet

    noct_max = new_df["nOctet"].max()
    noct_comp = new_df["nOctet"] == noct_max
    if noct_comp.sum() == 0:
        return "all failed"
    if noct_comp.sum() > 1:
        chgsepa_min = new_df[noct_comp]["ChgSepa"].min()
        chgsepa_comp = new_df[noct_comp]["ChgSepa"] == chgsepa_min
        if chgsepa_comp.sum() > 1:
            naro_max = new_df[noct_comp][chgsepa_comp]["nAromatic"].max()
            naro_comp = new_df[noct_comp][chgsepa_comp]["nAromatic"] == naro_max
            if naro_comp.sum() > 1:
                return ",".join(
                    new_df[noct_comp][chgsepa_comp][naro_comp].index.to_list()
                )
            else:
                return new_df[noct_comp][chgsepa_comp][naro_comp].index[0]
        else:
            return new_df[noct_comp][chgsepa_comp].index[0]
    else:
        return new_df[noct_comp].index[0]


cols = [
    "name",
    "chg",
    "ChgConserve",
    "NoRadical",
    "Sanitize",
    "BOSum",
    "ChgSepa",
    "nAromatic",
    "nOctet",
    "time",
]

obbLOG = "obabel/LOG.txt"
igxLOG = "indigox/LOG.txt"
x2mLOG = "xyz2mol/LOG.txt"
bosLOG = "bosolver/LOG.txt"
labLOG = "../label/LOG.txt"
# labLOG = "../label/LOG.proc"

# header = ["name","chg","ChgConserve","NoRadical","Sanitize","BOSum","ChgSepa","nAromatic","nOctet","time"]
obb = pd.read_csv(obbLOG, sep="\t", header=None)
obb.columns = cols

igx = pd.read_csv(igxLOG, sep="\t", header=None)
igx.columns = cols

x2m = pd.read_csv(x2mLOG, sep="\t", header=None)
x2m.columns = cols

bos = pd.read_csv(bosLOG, sep="\t", header=None)
bos.columns = cols

lab = pd.read_csv(labLOG, sep="\t", header=None)
lab.columns = cols

value_cols = [col for col in cols if col not in ["name", "chg"]]

for col in value_cols:
    obb[col] = obb[col].replace(-1, np.nan)
    igx[col] = igx[col].replace(-1, np.nan)
    x2m[col] = x2m[col].replace(-1, np.nan)
    bos[col] = bos[col].replace(-1, np.nan)
    lab[col] = lab[col].replace(-1, np.nan)  # ??
# print(bos[pd.isna(bos["time"])])

# print(f"{'Converted':=^20}")
# print(
#    f"OBB: {obb['time'].count()}\t "
#    + f"IGX: {igx['time'].count()}\t "
#    + f"X2M: {x2m['time'].count()}\t "
#    + f"BOS: {bos['time'].count()}"
# )
# print()
#
# print(f"LABEL: {len(lab)}")
# print()

obb = append_ChemVal_ChemStab(obb)
igx = append_ChemVal_ChemStab(igx)
x2m = append_ChemVal_ChemStab(x2m)
bos = append_ChemVal_ChemStab(bos)
lab = append_ChemVal_ChemStab(lab)

obb.to_pickle("df_obabel.pkl")
igx.to_pickle("df_indigox.pkl")
x2m.to_pickle("df_xyz2mol.pkl")
bos.to_pickle("df_bos.pkl")
# label does not need to be saved

assess_items = [
    "ChemVal",
    "ChgConserve",
    "NoRadical",
    "Sanitize",
    "ChemStab",
    "BOSum",
    "ChgSepa",
    "nAromatic",
    "nOctet",
]

df = pd.merge(
    obb[["name", "chg"] + assess_items],
    igx[["name", "chg"] + assess_items],
    on="name",
    suffixes=("_obb", "_igx"),
    how="outer",
)

df2 = pd.merge(
    x2m[["name", "chg"] + assess_items],
    bos[["name", "chg"] + assess_items],
    on="name",
    suffixes=("_x2m", "_bos"),
    how="outer",
)

df = pd.merge(
    df,
    df2,
    on="name",
    suffixes=("", ""),
    how="outer",
)

df = pd.merge(
    df,
    lab[["name", "chg"] + assess_items],
    on="name",
    how="outer",
)
df.rename(
    columns={
        "chg": "chg_lab",
        "ChgConserve": "ChgConserve_lab",
        "NoRadical": "NoRadical_lab",
        "Sanitize": "Sanitize_lab",
        "ChemVal": "ChemVal_lab",
        "ChemStab": "ChemStab_lab",
        "BOSum": "BOSum_lab",
        "ChgSepa": "ChgSepa_lab",
        "nAromatic": "nAromatic_lab",
        "nOctet": "nOctet_lab",
    },
    inplace=True,
)


c = df.columns.difference(["name"])
df[c] = df[c].astype("Int64")

df_val = df[["ChemVal_obb", "ChemVal_igx", "ChemVal_x2m", "ChemVal_bos", "ChemVal_lab"]]

val_obb = df_val["ChemVal_obb"] == 3
val_igx = df_val["ChemVal_igx"] == 3
val_x2m = df_val["ChemVal_x2m"] == 3
val_bos = df_val["ChemVal_bos"] == 3

print(f"{'ChemVal':=^20}")
print(
    f"OBB: {val_obb.sum()}\t IGX: {val_igx.sum()}\t X2M: {val_x2m.sum()}\t BOS: {val_bos.sum()}"
)
print()

val_lab = df_val["ChemVal_lab"] == 3
print(f"LABEL: {val_lab.sum()}")
print()

print(f"{'specific':-^20}")

ccs_obb = df["ChgConserve_obb"].sum()
ccs_igx = df["ChgConserve_igx"].sum()
ccs_x2m = df["ChgConserve_x2m"].sum()
ccs_bos = df["ChgConserve_bos"].sum()

nor_obb = df["NoRadical_obb"].sum()
nor_igx = df["NoRadical_igx"].sum()
nor_x2m = df["NoRadical_x2m"].sum()
nor_bos = df["NoRadical_bos"].sum()

san_obb = df["Sanitize_obb"].sum()
san_igx = df["Sanitize_igx"].sum()
san_x2m = df["Sanitize_x2m"].sum()
san_bos = df["Sanitize_bos"].sum()

print(
    pd.DataFrame(
        {
            "ChgConserve": [ccs_obb, ccs_igx, ccs_x2m, ccs_bos],
            "NoRadical": [nor_obb, nor_igx, nor_x2m, nor_bos],
            "Sanitize": [san_obb, san_igx, san_x2m, san_bos],
        },
        index=["OBB", "IGX", "X2M", "BOS"],
    )
)
print()
print("Unsuccessful (invalid)")
print(
    pd.DataFrame(
        {
            "ChgConserve": [
                df["ChgConserve_obb"].count() - ccs_obb,
                df["ChgConserve_igx"].count() - ccs_igx,
                df["ChgConserve_x2m"].count() - ccs_x2m,
                df["ChgConserve_bos"].count() - ccs_bos,
            ],
            "NoRadical": [
                df["NoRadical_obb"].count() - nor_obb,
                df["NoRadical_igx"].count() - nor_igx,
                df["NoRadical_x2m"].count() - nor_x2m,
                df["NoRadical_bos"].count() - nor_bos,
            ],
            "Sanitize": [
                df["Sanitize_obb"].count() - san_obb,
                df["Sanitize_igx"].count() - san_igx,
                df["Sanitize_x2m"].count() - san_x2m,
                df["Sanitize_bos"].count() - san_bos,
            ],
            "invalid": [
                (~val_obb).sum(),
                (~val_igx).sum(),
                (~val_x2m).sum(),
                (~val_bos).sum(),
            ],
        },
        index=["OBB", "IGX", "X2M", "BOS"],
    )
)
print()

#####################
#     CHEM STAB     #
#####################

# Set NaN to chemically invalid rows
stab_items = ["ChemStab", "BOSum", "ChgSepa", "nAromatic", "nOctet"]
for col in stab_items:
    # df[f"{col}_obb"][~val_obb] = np.nan
    df.loc[~val_obb, f"{col}_obb"] = np.nan
    # df[f"{col}_igx"][~val_igx] = np.nan
    df.loc[~val_igx, f"{col}_igx"] = np.nan
    # df[f"{col}_x2m"][~val_x2m] = np.nan
    df.loc[~val_x2m, f"{col}_x2m"] = np.nan
    # df[f"{col}_bos"][~val_bos] = np.nan
    df.loc[~val_bos, f"{col}_bos"] = np.nan
    # df[f"{col}_lab"][~val_lab] = np.nan
    df.loc[~val_lab, f"{col}_lab"] = np.nan

# save df as pickle
df.to_pickle("df_ana.pkl")


chemstab_res = df.apply(get_chemstab, axis=1)
chemstab_res["name"] = df["name"]
chemstab_res.to_pickle("chemstab_res.pkl")
general = [
    chemstab_res.str.contains("OBB").sum(),
    chemstab_res.str.contains("IGX").sum(),
    chemstab_res.str.contains("X2M").sum(),
    chemstab_res.str.contains("BOS").sum(),
]

print(f"{'ChemStab':=^20}")
print(f"OBB: {general[0]}\t IGX: {general[1]}\t X2M: {general[2]}\t BOS: {general[3]}")
print()
print("Unsuccessful (valid - chemstab)")
print(
    f"OBB: {val_obb.sum() - general[0]}\t IGX: {val_igx.sum() - general[1]}\t"
    + f" X2M: {val_x2m.sum() - general[2]}\t BOS: {val_bos.sum() - general[3]}"
)


cols = ["_obb", "_igx", "_x2m", "_bos"]
df_min_ChgSepa = df[["ChgSepa" + col for col in cols]].min(axis=1)
df_max_nAromatic = df[["nAromatic" + col for col in cols]].max(axis=1)
df_max_nOctet = df[["nOctet" + col for col in cols]].max(axis=1)

# NaN is not equal to NaN, so it is False in such case

# sta = [(df["ChemStab" + col] >= df_max) for col in cols]
chgsepa = [(df["ChgSepa" + col] <= df_min_ChgSepa) for col in cols]
naro = [(df["nAromatic" + col] >= df_max_nAromatic) for col in cols]
noct = [(df["nOctet" + col] >= df_max_nOctet) for col in cols]
# bochgarooct = [
#    (
#        (df["BOSum" + col] >= df_max_BOSum)
#        & (df["ChgSepa" + col] <= df_min_ChgSepa)
#        & (df["nAromatic" + col] >= df_max_nAromatic)
#        & (df["nOctet" + col] >= df_max_nOctet)
#    )
#    for col in cols
# ]


summary = pd.DataFrame(
    {
        #        "BOSum": [x.sum() for x in bosum],
        "nOctet": [x.sum() for x in noct],
        "ChgSepa": [x.sum() for x in chgsepa],
        "nAromatic": [x.sum() for x in naro],
    },
    # index=["OBB", "IGX", "X2M", "BOS", "LAB"],
    index=["OBB", "IGX", "X2M", "BOS"],
)

# print(f"OBB: {sta[0]}\t IGX: {sta[1]}\t X2M: {sta[2]}\t BOS: {sta[3]}")

# print()
# print(f"{'specific':-^20}")
print(summary)
print()

print("Unsuccessful (valid - chemstab)")
print(
    f"OBB: {val_obb.sum() - general[0]}\t IGX: {val_igx.sum() - general[1]}\t"
    + f" X2M: {val_x2m.sum() - general[2]}\t BOS: {val_bos.sum() - general[3]}"
)
# print summary of unsuccessful(valid - best)
# bosum_neg = [(df["BOSum" + col] < df_max_BOSum) for col in cols]
chgsepa_neg = [(df["ChgSepa" + col] > df_min_ChgSepa) for col in cols]
naro_neg = [(df["nAromatic" + col] < df_max_nAromatic) for col in cols]
noct_neg = [(df["nOctet" + col] < df_max_nOctet) for col in cols]
# bochgarooct_neg = [
#    (
#        (df["BOSum" + col] < df_max_BOSum)
#        | (df["ChgSepa" + col] > df_min_ChgSepa)
#        | (df["nAromatic" + col] < df_max_nAromatic)
#        | (df["nOctet" + col] < df_max_nOctet)
#    )
#    for col in cols
# ]

summary_neg = pd.DataFrame(
    {
        # "BOSum": [x.sum() for x in bosum_neg],
        "nOctet": [x.sum() for x in noct_neg],
        "ChgSepa": [x.sum() for x in chgsepa_neg],
        "nAromatic": [x.sum() for x in naro_neg],
        # "allBest": [x.sum() for x in bochgarooct_neg],
    },
    # index=["OBB", "IGX", "X2M", "BOS", "LAB"],
    index=["OBB", "IGX", "X2M", "BOS"],
)

print()
print(f"{'specific':-^20}")
print(summary_neg)
print()


# bosum_neg = [
#    ((df["BOSum" + col].isna()) | (df["BOSum" + col] < df_max_BOSum)) for col in cols
# ]
# chgsepa_neg = [
#    ((df["ChgSepa" + col].isna()) | (df["ChgSepa" + col] > df_min_ChgSepa))
#    for col in cols
# ]
# naro_neg = [
#    ((df["nAromatic" + col].isna()) | (df["nAromatic" + col] < df_max_nAromatic))
#    for col in cols
# ]
# noct_neg = [
#    ((df["nOctet" + col].isna()) | (df["nOctet" + col] < df_max_nOctet)) for col in cols
# ]
# bochgarooct_neg_na = [
#    (
#        (df["BOSum" + col].isna())
#        | (df["BOSum" + col] < df_max_BOSum)
#        | (df["ChgSepa" + col].isna())
#        | (df["ChgSepa" + col] > df_min_ChgSepa)
#        | (df["nAromatic" + col].isna())
#        | (df["nAromatic" + col] < df_max_nAromatic)
#        | (df["nOctet" + col].isna())
#        | (df["nOctet" + col] < df_max_nOctet)
#    )
#    for col in cols
# ]


# summary_neg_na = pd.DataFrame(
#    # {"allBest": [x.sum() for x in bochgarooct_neg_na]},
#    {"general": [x.sum() for x in bochgarooct_neg_na]},
#    # index=["OBB", "IGX", "X2M", "BOS", "LAB"],
#    index=["OBB", "IGX", "X2M", "BOS"],
# )
#
# print()
# print("Unsuccessful (total - best)")
# print(summary_neg_na)
# print()


t_obb_mean = obb["time"].dropna().mean()
t_obb_std = obb["time"].dropna().std()

t_igx_mean = igx["time"].dropna().mean()
t_igx_std = igx["time"].dropna().std()

t_x2m_mean = x2m["time"].dropna().mean()
t_x2m_std = x2m["time"].dropna().std()

t_bos_mean = bos["time"].dropna().mean()
t_bos_std = bos["time"].dropna().std()


print(f"{'TIME':=^20}")
print("OBB: {:.5f} ± {:.5f}".format(t_obb_mean, t_obb_std))
print("IGX: {:.5f} ± {:.5f}".format(t_igx_mean, t_igx_std))
print("X2M: {:.5f} ± {:.5f}".format(t_x2m_mean, t_x2m_std))
print("BOS: {:.5f} ± {:.5f}".format(t_bos_mean, t_bos_std))

print()
print(f"{'LONGER than 10s':-^20}")
print("OBB", (obb["time"] > 10).sum())
print("IGX", (igx["time"] > 10).sum())
print("X2M", (x2m["time"] > 10).sum())
print("BOS", (bos["time"] > 10).sum())

print()
df_t = pd.concat([obb["time"], igx["time"], x2m["time"], bos["time"]], axis=1)
df_t.columns = ["OBB", "IGX", "X2M", "BOS"]
# ax = sns.kdeplot(data=df_t, log_scale=True, fill=True)
colors = sns.color_palette().as_hex()
fig, axes = plt.subplots(4, 1, figsize=(10, 10))
axes[0].set_xlim((1e-5, 10))
axes[1].set_xlim((1e-5, 10))
axes[2].set_xlim((1e-5, 10))
axes[3].set_xlim((1e-5, 10))
sns.histplot(
    data=df_t[["OBB"]],
    x="OBB",
    log_scale=True,
    ax=axes[0],
    kde=False,
    binwidth=1e-1,
    color=colors[0],
)
sns.histplot(
    data=df_t[["IGX"]],
    x="IGX",
    log_scale=True,
    ax=axes[1],
    kde=False,
    binwidth=1e-1,
    color=colors[1],
)
sns.histplot(
    data=df_t[["X2M"]],
    x="X2M",
    log_scale=True,
    ax=axes[2],
    kde=False,
    binwidth=1e-1,
    color=colors[2],
)
sns.histplot(
    data=df_t[["BOS"]],
    x="BOS",
    log_scale=True,
    ax=axes[3],
    kde=False,
    binwidth=1e-1,
    color=colors[3],
)
axes[0].legend(["OBB"])
axes[1].legend(["IGX"])
axes[2].legend(["X2M"])
axes[3].legend(["BOS"])

for ax in axes:
    ax.set_ylabel(None)
fig.text(0.04, 0.5, "Count", va="center", rotation="vertical", fontsize=12)
plt.xlabel("Time (s)", fontsize=12)

plt.savefig("time.png")
#plt.show()
print("Time Distribution Plot was saved as time.png")

# print("OBB")
# print(df[~df["obb"]][["name", "ChemStab_obb", "max", "obb"]])
# print(df[df["ChemStab_obb"] - df_max < -5][["name", "ChemStab_obb", "max", "obb"]])

# print("IGX")
# print(df[~df["igx"]][["name", "ChemStab_igx", "max", "igx"]])
# print(df[df["ChemStab_igx"] - df["max"] < -5][["name", "ChemStab_igx", "max", "igx"]])
#
# print("X2M")
# print(df[~df["x2m"]][["name", "ChemStab_x2m", "max", "x2m"]])
# print(df[df["ChemStab_x2m"] - df["max"] < -5][["name", "ChemStab_x2m", "max", "x2m"]])
#
# print("BOS")
# print(df[~df["bos"]][["name", "ChemStab_bos", "max", "bos"]])
# print(df[df["ChemStab_bos"] - df["max"] < -5][["name", "ChemStab_bos", "max", "bos"]])

# compare(obb, igx, x2m, bos)
