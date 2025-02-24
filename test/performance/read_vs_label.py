import pandas as pd

df = pd.read_csv("vs_label.csv")
df.drop(columns=["Unnamed: 0"], inplace=True)

obb = df["obb"].value_counts()
igx = df["igx"].value_counts()
x2m = df["x2m"].value_counts()
ilp = df["ilp"].value_counts()

print("OBB")
print(obb)
print()
print("IGX")
print(igx)
print()
print("X2M")
print(x2m)
print()
print("ILP")
print(ilp)
print()
