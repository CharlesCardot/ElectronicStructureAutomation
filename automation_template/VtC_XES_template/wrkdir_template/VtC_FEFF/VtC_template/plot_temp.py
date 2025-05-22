import sys
import os
from pathlib import Path

path_to_file = Path(os.path.realpath(__file__))
automation_index = next((i for i, p in enumerate(path_to_file.parts) if p == 'automation'), None)
parent_dir = path_to_file.parents[len(path_to_file.parts) - automation_index - 3]

utils_path = parent_dir.parents[0] / "utils"
sys.path.append(str(utils_path))
import CharlesFunctions as CF

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("isotropic/xmu.dat", comments = "#").T
data = np.asarray([data[0],data[3]])
plt.plot(data[0],data[1],label="Iso")

data = np.loadtxt("x_polarization/xmu.dat", comments = "#").T
data = np.asarray([data[0],data[3]])
plt.plot(data[0],data[1],label="x pol")

data = np.loadtxt("y_polarization/xmu.dat", comments = "#").T
data = np.asarray([data[0],data[3]])
plt.plot(data[0],data[1],label="y pol")

data = np.loadtxt("z_polarization/xmu.dat", comments = "#").T
data = np.asarray([data[0],data[3]])
plt.plot(data[0],data[1],label="z pol")

plt.legend()
plt.title("NAME VtC-XES FEFF")
plt.show()




