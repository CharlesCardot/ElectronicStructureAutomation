import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
import re

def PlotTrim(arr, left, right):
    temp = [[],[]]
    for i in range(len(arr[0])):
        if left <= arr[0][i] <= right:
            temp[0].append(arr[0][i])
            temp[1].append(arr[1][i])
    return np.asarray(temp)

def readRIXSfile(filename):
    doc = list(open(filename, "r"))
    holder = []
    for i in range(len(doc)):
        if(doc[i] and doc[i][0] != '#' and i>4):
            txt = doc[i].split()
            for i in range(len(txt)):
                txt[i] = float(txt[i])
            holder.append(txt)
    x = np.transpose(np.asarray(holder))[0]
    temp = []
    yrange = int((len(holder[0])-3)/2)
    for i in range(yrange+1):
        temp2 = []
        for j in range(len(holder)):
            temp2.append(-1*holder[j][2+i*2])
        temp.append(np.asarray([x,temp2]))
    return np.asarray(temp)

def Convert_CompName_ToLaTeX(name):
    
    num_indexes = [m.start(0) for m in re.finditer('[0-9]+', name)]
    offset = 0
    for ind in num_indexes:
        name = name[0:ind+offset] + "$_\mathdefault{" + name[ind+offset] + "}$" + name[ind+offset+1:]
        name = r"{}".format(name)
        added_character_length = 17 # from '$_\mathdefault{}'
        offset += added_character_length
    return name

def makeplot_SingleRIXS(RIXS_filepath, omega1_min, omega1_max, plot_type, name, transpose=False, plot_params=None):
    '''

    ######################################################
    Incident => Excitation, for clarity (no photoelectron)
    ######################################################

    For plotting a Single RIXS.dat file. This relies on the output format from
    CreateResonantSpectra() in Quanty.
    '''
    
    plt.figure(figsize=(8,6))
    colormap = plt.get_cmap('magma_r')

    data = readRIXSfile(RIXS_filepath)
    spectra = np.asarray([i[1] for i in data])

    # x-axis
    x = data[0][0] 

    # y-axis
    f = open(RIXS_filepath,"r").readlines()
    NE1 = int(f[0].split()[-1])
    Emin1 = omega1_min
    Emax1 = omega1_max
    step = np.round((Emax1-Emin1)/(NE1-1),6)
    y = [Emin1 + step*i for i in range(len(spectra))]

    # z-matrix
    z=[]
    for i in range(len(spectra)):
        arr = np.asarray([x, spectra[i]])
        z.append(arr[1])
    
    plt.xlabel(r"Energy Loss, $\mathdefault{\omega_1} - \mathdefault{\omega_2}$ (eV)", fontsize=25)
    plt.ylabel(r"Incident Energy, $\mathdefault{\omega_1}$ (eV)", fontsize=25, labelpad=3)
    
    if transpose:
        x,y = y,x
        z = np.transpose(z)
        plt.xlabel(r"Incident Energy, $\mathdefault{\omega_1}}$ (eV)",fontsize=25)
        plt.ylabel(r"Energy Loss, $\mathdefault{\omega_1} - \mathdefault{\omega_2}$ (eV)", fontsize=25, labelpad=3)

    # Normalize peak height to 1 for nice colorbar
    z = z/np.max(z)

    # Green to yellow to orange to red
    colors = [(143 / 255, 197 / 255, 135 / 255), 
              (254 / 255, 232 / 255, 0 / 255), 
              (241 / 255, 159 / 255, 2 / 255), 
              (153 / 255, 36 / 255, 19 / 255)]
    color_positions = [0, 0.25, 0.60, 1.0]  # Corresponding to the colors defined
    custom_cmap = LinearSegmentedColormap.from_list("green_yellow_orange_red", list(zip(color_positions, colors)))

    # Draw contour plot, trace with black lines, add colorbar
    levels = np.linspace(0,1,20)
    contour = plt.contourf(x,y,z,levels,cmap=custom_cmap)
    contour_lines = plt.contour(x,y,z,levels,colors='black',linewidths=0.5)
    plt.colorbar(contour)
    
    ## Old
    #cbar = plt.colorbar(im)
    #cbar = plt.colorbar(im, ticks=[0.0,0.25,0.5,0.75,1.0])
    #tick_font_size = 18
    #cbar.ax.tick_params(labelsize=tick_font_size)
        
    text = Convert_CompName_ToLaTeX(name) + " RIXS"
    plt.gca().text(0.05, 0.95, text, ha='left', va='top', fontsize = 25, transform=plt.gca().transAxes)
    plt.show()


def makeplot_SingleRIXS_fullcuts(RIXS_filepath, omega1_min, omega1_max, plot_type, name, transpose=False, plot_params=None):
    '''

    ######################################################
    Incident => Excitation, for clarity (no photoelectron)
    ######################################################

    For plotting a Single RIXS.dat file. This relies on the output format from
    CreateResonantSpectra() in Quanty.
    '''
    
    fig, axs = plt.subplots(2, 2, figsize=(8,6))

    colormap = plt.get_cmap('magma_r')

    data = readRIXSfile(RIXS_filepath)
    spectra = np.asarray([i[1] for i in data])

    # x-axis (Energy Loss)
    x = data[0][0] 

    # y-axis (Incident Energy)
    f = open(RIXS_filepath,"r").readlines()
    NE1 = int(f[0].split()[-1])
    Emin1 = omega1_min
    Emax1 = omega1_max
    step = np.round((Emax1-Emin1)/(NE1-1),6)
    y = [Emin1 + step*i for i in range(len(spectra))]

    #############################################
    # Transpose by default 
    # Incident goes to Energy Loss and vice versa
    #############################################
    x,y = y,x
    spectra = np.transpose(spectra)
    axs[1][0].set_xlabel(r"Incident Energy, $\mathdefault{\omega_1}}$ (eV)", fontsize=15)
    axs[1][0].set_ylabel(r"Energy Loss, $\mathdefault{\omega_1} - \mathdefault{\omega_2}$ (eV)", fontsize=15, labelpad=3)

    # Normalize peak height to 1 for nice colorbar
    z = spectra / np.max(spectra)

    # Green to yellow to orange to red
    colors = [(143 / 255, 197 / 255, 135 / 255), 
              (254 / 255, 232 / 255, 0 / 255), 
              (241 / 255, 159 / 255, 2 / 255), 
              (153 / 255, 36 / 255, 19 / 255)]
    color_positions = [0, 0.25, 0.60, 1.0]  # Corresponding to the colors defined
    custom_cmap = LinearSegmentedColormap.from_list("green_yellow_orange_red", list(zip(color_positions, colors)))

    # Draw contour plot, trace with black lines, add colorbar
    levels = np.linspace(0,1,20)
    contour = axs[1][0].contourf(x, y, z, levels, cmap=custom_cmap)
    contour_lines = axs[1][0].contour(x,y,z,levels,colors='black',linewidths=0.5)
    cbar = plt.colorbar(contour, ax=axs[1][0])
    cbar.set_ticks([np.min(z), np.max(z)])
    cbar.set_ticklabels(['min', 'max'])
    
    text = Convert_CompName_ToLaTeX(name) + " RIXS"
    fig.suptitle(text, y=0.96, fontsize = 20)

    ##########################
    # Constant Energy Transfer
    ##########################

    y_cet = np.sum(z, axis=0)
    axs[0][0].plot(x, y_cet)
    axs[0][0].text(0.90, 0.90, "CET", transform = axs[0][0].transAxes, va='top', ha='right', fontsize = 15)


    ##########################
    # Constant Incident Energy
    ##########################

    y_cie = np.sum(z, axis=1)
    axs[1][1].plot(y, y_cie)
    axs[1][1].text(0.90, 0.90, "CIE", transform = axs[1][1].transAxes, va='top', ha='right', fontsize = 15)

    ##########################
    # Constant Emission Energy
    ##########################

    spectra_points = []
    """
    x = Incident Energy (2000)
    y = Energy Transfer (201)
    
    """
    for j in range(len(spectra)):
        for i in range(len(spectra[0])):
            spectra_points.append((x[i], (y[j] - x[i]) * -1, spectra[j][i]))

    spectra_x = [i[0] for i in spectra_points]
    spectra_y = [i[1] for i in spectra_points]
    spectra_z = [i[2] for i in spectra_points]

    # Create a meshgrid from x and y
    x_mesh = np.linspace(min(spectra_x), max(spectra_x), 100)
    y_mesh = np.linspace(min(spectra_y), max(spectra_y), 100)
    X, Y = np.meshgrid(x_mesh, y_mesh)

    # Interpolate z values on the meshgrid points
    Z = np.zeros_like(X)
    for point_x, point_y, point_z in spectra_points:
        Z += np.exp(-((X - point_x) ** 2 + (Y - point_y) ** 2) / (2 * 0.1 ** 2)) * point_z
    Z = Z / np.max(Z)

    # Create a contour plot
    y_cee = np.sum(Z, axis=0)
    axs[0][1].plot(x_mesh, y_cee)
    
    axs[0][1].text(0.90, 0.90, "CEE", transform = axs[0][1].transAxes, va='top', ha='right', fontsize = 15)

    plt.show()



import json
with open("spectra_params.json", "r") as f:
    spectra_params = json.load(f)

RIXS_filepath = "RIXS.dat"
omega1_min = spectra_params["Emin1"]
omega1_max = spectra_params["Emax1"]
plot_type = "contour"
name = "MATERIAL_NAME"

makeplot_SingleRIXS_fullcuts(RIXS_filepath, omega1_min, omega1_max, plot_type, name)
