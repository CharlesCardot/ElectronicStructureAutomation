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

data = np.loadtxt("XES_xpol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
x_shift = data[0][np.argmax(data[1])]
data[0] = data[0] - x_shift
data = CF.PlotTrim(CF.Flipping(data),left,right)
plt.plot(data[0],data[1],label="x pol")

data = np.loadtxt("XES_ypol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data[0] = data[0] - x_shift
data = CF.PlotTrim(CF.Flipping(data),left,right)
plt.plot(data[0],data[1],label="y pol")

data = np.loadtxt("XES_zpol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data[0] = data[0] - x_shift
data = CF.PlotTrim(CF.Flipping(data),left,right)
plt.plot(data[0],data[1],label="z pol")

plt.legend()
plt.title("NAME 3p->1s XES Quanty")
plt.show()




