import sys
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF

from scipy.interpolate import CubicSpline
from scipy import optimize
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os

left = -10
right = 25

data = np.loadtxt("XAS_xpol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data = CF.PlotTrim(data,left,right)
data = CF.PeakNormalize(data)
plt.plot(data[0],data[1],label="x_pol 2p XAS Quanty")

data = np.loadtxt("XAS_ypol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data = CF.PlotTrim(data,left,right)
data = CF.PeakNormalize(data)
plt.plot(data[0],data[1],label="y_pol 2p XAS Quanty")

data = np.loadtxt("XAS_zpol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data = CF.PlotTrim(data,left,right)
data = CF.PeakNormalize(data)
plt.plot(data[0],data[1],label="z_pol 2p XAS Quanty")

plt.xlabel('Energy (eV)', fontsize = 20)
plt.ylabel('Intensity (arb. units)', fontsize = 20)

plt.tight_layout()
plt.legend()
plt.show()




