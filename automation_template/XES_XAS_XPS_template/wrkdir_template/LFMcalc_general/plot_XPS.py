import sys
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF
import numpy as np
import matplotlib.pyplot as plt
import os

left = -50
right = 50

data = np.loadtxt("XPS.dat",skiprows=5).T
data = np.asarray([data[0],-1*data[2]])
data = CF.PlotTrim(data,left,right)
data = CF.PeakNormalize(data)
plt.plot(data[0],data[1])

plt.xlabel('Energy (eV)', fontsize = 20)
plt.ylabel('Intensity (arb. units)', fontsize = 20)
plt.title('NAME 2p XPS', fontsize = 20)

plt.tight_layout()
plt.legend()
plt.show()




