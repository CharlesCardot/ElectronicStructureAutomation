import sys
import traceback
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os
import re
import time

from scipy.signal import find_peaks

# Regular expression to match the folder format D_x+U_x
pattern = re.compile(r'D_([-]?\d{1,3}\.\d{1})\+U_([-]?\d{1,3}\.\d{1})')

# Set to store unique D and U values
Delta, Udd = set(), set()

# Search through the current directory
for folder in os.listdir('.'):
    match = pattern.match(folder)
    if match:
        D, U = float(match.group(1)), float(match.group(2))
        Delta.add(D)
        Udd.add(U)

# Sort the unique values, sorted(set) returns a list
Delta = sorted(Delta)
Udd = sorted(Udd)

left = -5
right = 20

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(13.5,4.5))

Delta = [2.0]
# Udd = [5.0]
for D in range(len(Delta)):
    for U in range(len(Udd)):
        XES_filename = "D_"+str(Delta[D])+"+U_"+str(Udd[U])+"/XES.dat"
        XAS_filename = "D_"+str(Delta[D])+"+U_"+str(Udd[U])+"/XAS.dat"
        XPS_filename = "D_"+str(Delta[D])+"+U_"+str(Udd[U])+"/XPS.dat"
        
        try:
            # PLOTTING BROADENED XES
            data = np.loadtxt(XES_filename,skiprows=5).T
            data = np.asarray([data[0],-1*data[2]])
            data = CF.alignpeak(data)
            data = CF.Flipping(data)
            data = CF.PlotTrim(data, -25, 5)
            
            if len(Udd) < 2:
                label = r'$\Delta$ = ' + f'{Delta[D]}'
                legend_title = r"U$_{dd}$ = " + str(Udd[U])
            elif len(Delta) < 2: 
                label = r'U$_{dd}$ = ' + f'{Udd[U]}'
                legend_title = r"$\Delta$ = " + str(Delta[D])
            else:
                legend_title = ''
                label = r'$\Delta$ = ' + f'{Delta[D]}, ' + r'U$_{dd}$ = ' + f'{Udd[U]}'

            ax[0].plot(data[0], data[1], label=label)
            ax[0].set_title('XES', fontsize = 20)
 
            # PLOTTING BROADENED XAS
            data = np.loadtxt(XAS_filename,skiprows=5).T
            data = np.asarray([data[0],-1*data[2]])
            data = CF.alignpeak(data)
            data = CF.PlotTrim(data, -5, 35)
            
            ax[1].plot(data[0], data[1], label=label)
            ax[1].set_title('XAS', fontsize = 20)
            
            # PLOTTING BROADENED XPS
            data = np.loadtxt(XPS_filename,skiprows=5).T
            data = np.asarray([data[0],-1*data[2]])
            data = CF.alignpeak(data)
            data = CF.PlotTrim(data, -5, 35)
            
            ax[2].plot(data[0], data[1], label=label)
            ax[2].set_title('XPS', fontsize = 20)
 
        except Exception:
            print("Could not load", "D_"+str(Delta[D])+"+U_"+str(Udd[U]))
            print(traceback.format_exc())
        
    ax[0].legend(fontsize=10, title=legend_title)
    ax[0].set_ylabel("Intensity (arb. units)", fontsize=15)
    fig.suptitle("CuO XAS, XPS, and XES", fontsize=20, y=0.98)
    plt.tight_layout()
    # plt.savefig("NAME_CT_param_explore.png", dpi=100)
    # print("Successfully created", "CuO_CT_param_explore.png")
    plt.show()





























