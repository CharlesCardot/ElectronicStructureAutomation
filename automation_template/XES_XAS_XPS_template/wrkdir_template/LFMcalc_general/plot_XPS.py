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
plt.show()




