import sys
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF
import numpy as np
import matplotlib.pyplot as plt
import os

left = -25
right = 10

data = np.loadtxt("XES_xpol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data = CF.PlotTrim(CF.Flipping(CF.alignpeak(data)),left,right)
#data = CF.PeakNormalize(data)
plt.plot(data[0],data[1],label="x_pol")

data = np.loadtxt("XES_ypol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data = CF.PlotTrim(CF.Flipping(CF.alignpeak(data)),left,right)
#data = CF.PeakNormalize(data)
plt.plot(data[0],data[1],label="y_pol")

data = np.loadtxt("XES_zpol.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data = CF.PlotTrim(CF.Flipping(CF.alignpeak(data)),left,right)
#data = CF.PeakNormalize(data)
plt.plot(data[0],data[1],label="z_pol")

plt.xlabel('Energy (eV)', fontsize = 20)
plt.ylabel('Intensity (arb. units)', fontsize = 20)
plt.title('NAME Kalpha XES', fontsize = 20)

plt.tight_layout()
plt.legend()
plt.show()




