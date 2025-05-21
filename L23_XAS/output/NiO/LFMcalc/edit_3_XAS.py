import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from pathlib import Path

path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]

utils_path = parent_dir.parents[0] / "utils"
sys.path.append(str(utils_path))
import CharlesFunctions as CF



with open("testing.out","r") as f:
    data = [x.strip() for x in f.readlines()]

Energies = CF.badGREP(data,start_str="--Initial State Energies--",start=1,end=2)[0]
Energies = [float(x.replace("{","").replace("}","").strip()) for x in Energies.split(",")]
Energies = np.asarray(Energies)

Egnd = Energies[0]
Npsi = 1
T = 293/11606
for Energy in Energies:
    dZ = np.exp(-(Energy - Egnd)/T)
    Npsi += 1
    if dZ < 0.001:
        break
print("Npsi", Npsi)


rho_d = CF.badGREP(data,start_str="Calculate the DFT 1-particle density matrix",start=2,end=2+10)
rho_d[0] = "rho_d = " + rho_d[0]
rho_d = [x + "\n" for x in rho_d]
print("rho_d", rho_d)

# Edit 3_XAS.Quanty to make it run faster
filename = "3_XAS.Quanty"
with open(filename, "r") as f:
    lines = f.readlines()
    
rho_line_offset = None
for key, line in enumerate(lines):

    if "Npsi=" in line:
        line = "Npsi=" + str(Npsi) + "\n"
        lines[key] = line
    if "Next we need to correct for the double counting" in line:
        rho_line_offset = key + 10
        for i in range(key, key + 10):
            lines[i] == "--" + lines[i]

lines = lines[0:rho_line_offset] + rho_d + lines[rho_line_offset:]

with open(filename,"w") as f:
    for line in lines:
        f.write(line)






