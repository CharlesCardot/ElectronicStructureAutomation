import sys
import traceback
sys.path.insert(1,"/home/ccardot3/Python_Code/CharlesFunctions/")
import CharlesFunctions as CF

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy
import os
import time

from scipy.signal import find_peaks

Udd = [2.0, 4.0, 6.0, 8.0]
Delta = ["{:.1f}".format(x) for x in np.arange(-10,10.1,2)]

left = -5
right = 20

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(13.5,4.5), sharex=True)
fig.subplots_adjust(wspace=0)

Udd = [2.0]
for U in range(len(Udd)):
    x_normal = 0
    x_sticks = 0
    for D in range(len(Delta)):
        name_normal = "D_"+str(Delta[D])+"+U_"+str(Udd[U])+"/XES.dat"
        name_sticks = "D_"+str(Delta[D])+"+U_"+str(Udd[U])+"/XES_sticks.dat"
        
        try:
            # PLOTTING BROADENED
            data = np.loadtxt(name_normal,skiprows=5).T
            data = np.asarray([data[0],-1*data[2]])
            data = CF.PeakNormalize(CF.PlotTrim(CF.alignpeak(data),left,right))
            data = CF.Flipping(data)
            x_normal = x_normal - np.max(data[1])*0.1
            
            # FOR FWHM
            #peakpos = find_peaks(data[1],prominence=(0.01,1))[0]
            #FWHM = np.round(CF.FWHM(data,peakpos[1]),4)
            #label = r"$\Delta = $" + str(Delta[D]) + ", FWHM: " + str(FWHM)

            label = r"$\Delta = $" + str(Delta[D])

            ax[0].plot(data[0],data[1],label=label)
            ax[1].plot(data[0],data[1]+x_normal,label=label)
            ax[1].set_yticks([])

            # PLOTTING STICKS
            data = np.loadtxt(name_sticks,skiprows=5).T
            data = np.asarray([data[0],-1*data[2]])
            data = CF.PeakNormalize(CF.PlotTrim(CF.alignpeak(data),left,right))
            data = CF.Flipping(data)
            x_sticks = x_sticks - np.max(data[1])*0.1

            ax[2].plot(data[0],data[1]+x_sticks,label=label)
            ax[2].set_yticks([])

        except Exception:
            print("Could not load", "D_"+str(Delta[D])+"+U_"+str(Udd[U]))
            print(traceback.format_exc())
        
    legend_title = r"U$_{dd}$ " + str(Udd[U])
    ax[0].legend(fontsize=10, title=legend_title)
    ax[0].set_ylabel("Intensity (arb. units)", fontsize=15)
    fig.suptitle("NAME", fontsize=20, y=0.95)
    plt.savefig("NAME_CT_param_explore.png", dpi=100)
    print("Successfully created", "NAME_CT_param_explore.png")
    #plt.show()





























