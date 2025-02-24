import pandas as pd
import pickle


def compare(r1, r2):
    if r1["ChemVal"] is pd.NA and r2["ChemVal"] is pd.NA:
        return 0
    elif r1["ChemVal"] is pd.NA:
        return -1
    elif r2["ChemVal"] is pd.NA:
        return 1

    if r1["ChemVal"] != 3 and r2["ChemVal"] != 3:
        return 0
    if r1["ChemVal"] > r2["ChemVal"]:
        return 1
    elif r1["ChemVal"] < r2["ChemVal"]:
        return -1
    else:
        if r1["nOctet"] > r2["nOctet"]:
            return 1
        elif r1["nOctet"] < r2["nOctet"]:
            return -1
        else:
            if r1["ChgSepa"] < r2["ChgSepa"]:
                return 1
            elif r1["ChgSepa"] > r2["ChgSepa"]:
                return -1
            else:
                if r1["nAromatic"] > r2["nAromatic"]:
                    return 1
                elif r1["nAromatic"] < r2["nAromatic"]:
                    return -1
                else:
                    return 0


def compare_to_label(row: pd.Series):
    name = row["name"]
    cols = ["ChemVal", "nOctet", "ChgSepa", "nAromatic"]
    obb = row[[col + "_obb" for col in cols]]
    igx = row[[col + "_igx" for col in cols]]
    x2m = row[[col + "_x2m" for col in cols]]
    bos = row[[col + "_bos" for col in cols]]
    lab = row[[col + "_lab" for col in cols]]

    obb.index = cols
    igx.index = cols
    x2m.index = cols
    bos.index = cols
    lab.index = cols

    res = {
        "name": name,
        "obb": compare(obb, lab),
        "igx": compare(igx, lab),
        "x2m": compare(x2m, lab),
        "bos": compare(bos, lab),
    }
    return res


if __name__ == "__main__":
    with open("df_ana.pkl", "rb") as f:
        df_ana = pickle.load(f)

    compared = df_ana.apply(compare_to_label, axis=1, result_type="expand")
    compared.to_csv("vs_label.csv")
