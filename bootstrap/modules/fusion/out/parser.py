import re

filename = "candidates_B_DE.txt"
vals = []
r    = []
with open(filename, mode="r") as file:
    while True:
        line = file.readline()
        if not line:
            break

        if "  b^1 = " in line:
            vals.append(float(re.search(r'\[\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line).group(1)))

        if "  R =" in line:
            r.append(float(re.search(r'\[\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line).group(1)))

print(len(vals), sorted(set(vals)))
print(len(r),sorted(set(r)))
