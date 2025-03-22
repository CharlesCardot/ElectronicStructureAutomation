import os
import sys

sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF
from pathlib import Path

import numpy as np



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
# This will speed up the batch jobs
filename = "3_XAS.Quanty"
with open(filename, "r") as f:
    lines = f.readlines()
    
for key, line in enumerate(lines):

    if "Npsi=" in line:
        line = "Npsi=" + str(Npsi) + "\n"
    if 220 < key < 230:
        line = "--" + line
        
    lines[key] = line
lines = lines[0:230] + rho_d + lines[230:]

with open(filename,"w") as f:
    for line in lines:
        f.write(line)






