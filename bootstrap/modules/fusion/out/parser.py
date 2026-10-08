import re

filename = "candidates_B_DE.txt"
vals = []
with open(filename, mode="r") as file:
    while True:
        line = file.readline()
        if not line:
            break

        if "  |b^1| = " in line:
            vals.append(float(re.search(r'\[\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)', line).group(1)))

print(sorted(set(vals)))
