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

left = -25
right = 10

data = np.loadtxt("XAS.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data = CF.PlotTrim(CF.Flipping(CF.alignpeak(data)),left,right)
data = CF.PeakNormalize(data)
plt.plot(data[0],data[1],label="NAME 2p XAS Quanty")

plt.legend()
plt.show()




