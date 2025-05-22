import sys
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF

import numpy as np
import matplotlib.pyplot as plt
import scipy
import os

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




