import os
import sys

sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF
from pathlib import Path

import numpy as np



with open("testing.out","r") as f:
    data = [x.strip() for x in f.readlines()]

rho_d = CF.badGREP(data,start_str="Calculate the DFT 1-particle density matrix",start=2,end=2+10)
rho_d[0] = "rho_d = " + rho_d[0]
rho_d = [x + "\n" for x in rho_d]
print("rho_d", rho_d)

# Edit script.Quanty to make it run faster
filename = "script.Quanty"
with open(filename, "r") as f:
    lines = f.readlines()
    
rho_line_offset = None
for key, line in enumerate(lines):

    if "Next we need to correct for the double counting" in line:
        rho_line_offset = key + 10
        for i in range(key, key + 10):
            lines[i] == "--" + lines[i]

lines = lines[0:rho_line_offset] + rho_d + lines[rho_line_offset:]

with open(filename,"w") as f:
    for line in lines:
        f.write(line)

#################################################################################
# Change qsub.script so that it no longer contains the call to edit script.Quanty
filename = "qsub.script"
with open(filename, "r") as f:
    lines = f.readlines()

for key, line in enumerate(lines):
    lines[key] = line.replace('python', '# python')

with open(filename,"w") as f:
    for line in lines:
        f.write(line)




