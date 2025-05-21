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

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import scipy
import shutil
import datetime
import os

from PyPDF2 import PdfWriter
from PIL import Image, ImageDraw, ImageFont
from PIL import ImageOps

from pymatgen.vis.structure_vtk import StructureVis
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer



"""
Notation suggested by Seidler:
"Suggestion for notation;  (prop dir, linear pol dir (can be zero), circular pol helicity (can be zero)).
Hence a spectrum of the z-directed XES with x-polarization would be (z,x,0) 
while a spectrum for the x-directed right-circular XES would be (x,0,+1)."
"""


def extract_quanty_spectra(filename):
    data = np.loadtxt(filename, skiprows = 5).T
    data = np.asarray([data[0],-1*data[2]])
    return data

def extract_feff_spectra(filename):
    data = np.loadtxt(filename, comments = "#").T
    data = np.asarray([data[0], data[3]])
    return data

def process_quanty_spectra_xes(data, x_shift = 0, left = -100, right = 100):
    data[0] += x_shift
    data = CF.PlotTrim(CF.Flipping(data), left, right)
    return data

def process_quanty_spectra_xas(data, x_shift = 0, left = -100, right = 100):
    data[0] += x_shift
    data = CF.PlotTrim(data, left, right)
    return data

def process_feff_spectra(data, x_shift = 0, left = -100, right = 100):
    data[0] += x_shift
    data = CF.PlotTrim(data, left, right)
    return data

def avg_two_spectra(dataone, datatwo):
    if not np.all(np.abs(dataone[0] - datatwo[0]) < 1e-14):
        error_msg = """x-axis of the spectra do not match, 
        averaging has not been implement for this case"""
        raise ValueError(error_msg)
    data = np.asarray([dataone[0], 0.5 * (dataone[1] + datatwo[1])])
    return data

def make_SCFConvPlots(fplo_out_path, plot_params_dict):

    with open(Path(fplo_out_path), "r") as f:
        lines = f.readlines()
    
    scf_lines = [line for line in lines if "SCF: iteration" in line]
    scf_lines = scf_lines[1:] # The first SCF: iteration isn't actual data
    scf_energies = []
    converged = False
    for line in scf_lines:
        if "CONVERGED" in line:
            converged = True
        line = line.split("u=")[-1]
        energy = float(line.split()[0])
        scf_energies.append(energy)
     
    x = [i for i in range(len(scf_energies))]
    plt.plot(x, scf_energies)
    plt.scatter(x, scf_energies)

    if converged:
        plt.axhline(y=scf_energies[-1], color="gray", linewidth = 2, dashes=[2,2])
        plt.text(0, scf_energies[-1] * 2, "Converged", color = "tab:green", fontsize = 12, ha="left")
    else:
        plt.text(0, scf_energies[-1] * 2, "NOT CONVERGED", color = "tab:red", fontsize = 15, ha="left")
    
    plt.yscale("log")
    plt.grid(which='both', linestyle='-', linewidth=0.5)
    plt.gca().yaxis.grid(True, which='minor', linestyle='--')
    
    plt.xticks(x)
    plt.xlabel("Iteration Step", fontsize=20)
    plt.ylabel("Energy (Hartree)", fontsize=20)

    plt.title("SCF Convergence, " + str(Path(fplo_out_path).name), fontsize=20)
    plt.tight_layout()

    plt.savefig(plot_params_dict["output_filename"], dpi=150, bbox_inches='tight')
    plt.close()

def make_BandStructAndLdosPlots(dft_path, material_input, plot_params_dict):
    """
    Plot the band structure, TM 3d ldos, and Ligand valence ldos.
    """

    def readfile(filename):
        doc = list(open(filename, "r"))
        holder = []
        for i in range(len(doc)):
            txt = doc[i].split()
            holder.append(txt)
        #return np.asarray(holder)
        return holder
    
    def readband(filename):
        data = np.loadtxt(filename).T
        x = np.asarray(data[0])
    
        temp = []
        for i in range(len(data)):
            temp.append(data[i])
    
        temp = np.asarray(temp)
        return x, temp
    
    def readpoints(filename):
        doc = list(open(filename, "r"))
        holder = []
        for i in range(len(doc)):
            txt = doc[i].split()
            holder.append(txt)
        holder = [x for x in holder if x]
        temp = []
        for i in range(len(holder)):
            if(i != 0 and holder[i][0] == '#'):
                if(holder[i][2] == '$~G'):
                    temp.append(np.asarray(["G",float(holder[i+1][0])]))
                else:
                    temp.append(np.asarray([holder[i][2],float(holder[i+1][0])]))
        return np.asarray(temp)
    
    def readwanbweights(filename):
        headerlist = list(open(filename, "r"))[1].split("  ")
        headerlist = [x.strip() for x in headerlist if x]
        headerlist = [x for x in headerlist if x]
    
    
        data = np.loadtxt(filename)
        temp = data.T[0]
        bands = len(data.T[2:])
    
        res = []
        for i in temp:
            if res.count(i) == 0:
                res.append(i)
    
        x = np.asarray(res)
        y = np.reshape(data.T[1],(-1,bands)).T
    
        wanbweights = {}
        for i in range(bands):
            wanbweights[headerlist[i+3]] = np.reshape(data.T[2+i],(-1,bands)).T
    
        return x, y, wanbweights
    
    def readbweights(filename,TMname,Ligandname):
        headerlist = list(open(filename, "r"))[1].split("  ")
        headerlist = [x.strip().replace(" ","") for x in headerlist]
        headerlist = [x for x in headerlist if x]
    
        data = np.loadtxt(filename)
        temp = data.T[0]
        bands = len(data.T[2:])
    
        res = []
        for i in temp:
            if res.count(i) == 0:
                res.append(i)
    
        x = np.asarray(res)
        y = np.reshape(data.T[1],(-1,bands)).T
    
        bweights = {}
        for i in range(bands):
            bweights[headerlist[i+3]] = np.reshape(data.T[2+i],(-1,bands)).T
    
        return x, y, bweights
    
    def readwandos(dft_path,TM,TMorb,Ligand,Ligandorb):
    
        #Actually reads in the wan l DOS
        TMtemp = []
        Ligandtemp = []
        for filename in os.listdir(dft_path):
            if filename.startswith("+wanldos."):
                doc = readfile(dft_path / filename)
                header = ' '.join([str(elem) for elem in doc[0]])
                if(TMorb in header and TM in header):
                    TMtemp.append(np.loadtxt(dft_path / filename).T[1])
                elif(Ligandorb in header and Ligand in header):
                    Ligandtemp.append(np.loadtxt(dft_path / filename).T[1])
    
        TMdos = np.sum(TMtemp,axis=0)
        Liganddos = np.sum(Ligandtemp,axis=0)
    
    
        x = np.loadtxt(dft_path / "+wanldos.001.spin1").T[0]
    
        return x, TMdos, Liganddos
    
    
    def readdos(dft_path,TM,TMorb,Ligand,Ligandorb):
        doc = readfile(dft_path / "=.pipe")
        temp = {}
        for i in range(len(doc)):
            if doc[i] == list(['#', 'Now,', 'give', 'list', 'of', 'ALL', '!!!', 'Wyckoff', 'positions.']):
                while(doc[i+1][0] != '#'):
                    temp[doc[i+1][0][1]]=doc[i+1][1]
                    i = i+1
    
        TMindexes = [k for k,v in temp.items() if v == TM]
        Ligandindexes = [k for k,v in temp.items() if v == Ligand]
    
        #Actually reads in the DOS
        temp = []
        for i in TMindexes:
            for filename in os.listdir(dft_path):
                if filename.startswith("+dos.sort00"+str(i)+".nl"):
                    doc = readfile(dft_path / filename)
                    if(doc[0][6]==TMorb):
                        temp.append(np.loadtxt(dft_path / filename).T[1])
        TMdos = np.sum(temp,axis=0)
    
        temp = []
        for i in Ligandindexes:
            for filename in os.listdir(dft_path):
                if filename.startswith("+dos.sort00"+str(i)+".nl"):
                    doc = readfile(dft_path / filename)
                    if(doc[0][6]==Ligandorb):
                        temp.append(np.loadtxt(dft_path / filename).T[1])
        Liganddos = np.sum(temp,axis=0)
    
        x = np.loadtxt(dft_path / "+dos.sort001.nl001").T[0]
    
        return x, TMdos, Liganddos
    
    
    def orbitalbands(weights, TM, TMorb, Ligand, Ligandorb):
        """
        Orbital bands are basically just a single curve on the bandstructure plot.
        All the ik values for where in reciprocal space the band lives has been parsed
        and a single vale in TMweights for example is the total weight (percentage TM)
        that a specific band in k-space has.
        """

        weightskeys = list(weights.keys())
        TMkeys = []
        Ligandkeys = []
    
        for i in weightskeys:
            if TM in i and TMorb in i:
                TMkeys.append(i)
            elif Ligand in i and Ligandorb in i:
                Ligandkeys.append(i)
    
        TMweights = []
        Ligandweights = []
    
        for i in TMkeys:
            TMweights.append(weights[i])
        TMweights = np.sum(TMweights,axis=0)
    
        for i in Ligandkeys:
            Ligandweights.append(weights[i])
        Ligandweights = np.sum(Ligandweights,axis=0)
    
        TMtoplot = []
        for i in range(len(TMweights)):
            if np.any(TMweights[i]>0.5):
                TMtoplot.append(i)
    
        Ligandtoplot = []
        for i in range(len(Ligandweights)):
            if np.any(Ligandweights[i]>0.5):
                Ligandtoplot.append(i)
    
        return TMtoplot, Ligandtoplot 

    def makeBandStructDOSplot(dft_path, output_filename, emin, emax, compound_name, TM, Ligand, TMorb = "3d", Ligandorb = "2p"):
    
        if not(TM in compound_name and Ligand in compound_name):
            raise ValueError("TM and Ligand not in compound name, check your definitions!")
    
        xwan, ywan, wanbweights = readwanbweights(dft_path / "+wanbweights")
        x, y, _ = readbweights(dft_path / "+bweights", TM, Ligand)
        highsympoints = readpoints(dft_path / "+points")
    
        highsympointlabels = highsympoints.T[0]
        highsympoints = highsympoints.T[1]
        highsympoints = [float(x) for x in highsympoints]
    
        TMcolor = "tab:red"
        Ligandcolor = "tab:blue"
        TMtoplot, Ligandtoplot = orbitalbands(wanbweights,TM,TMorb,Ligand,Ligandorb)
    
        fig, ax = plt.subplots(1,2,figsize=(8,5),gridspec_kw={'width_ratios': [3, 1]})
        plt.subplots_adjust(wspace=0.01)
    
        ##Bandstructure plot
        #Plotting Bands
        for i in range(len(TMtoplot)):
            ax[0].plot(xwan,ywan[TMtoplot[i]],linewidth=3,color=TMcolor,alpha=0.8)
        for i in range(len(Ligandtoplot)):
            ax[0].plot(xwan,ywan[Ligandtoplot[i]],linewidth=3,color=Ligandcolor,alpha=0.8)
        for i in range(len(TMtoplot)):
           for j in range(len(Ligandtoplot)):
               if(TMtoplot[i]==Ligandtoplot[j]):
                   ax[0].plot(xwan,ywan[Ligandtoplot[j]],linewidth=3,color="tab:purple",alpha=1)
 
    
        #Plotting High Symmetry Points (ex: G, E, K)
        for i in y:
            ax[0].plot(x,i,color="black",linewidth=0.5)
        for i in highsympoints:
            ax[0].vlines(float(i),emin,emax,linewidth=0.5,color="black")
    
        ax[0].set_ylim((emin,emax))
        ax[0].set_xlim((x[0],x[-1]))
        ax[0].set_xticks(highsympoints)
        ax[0].set_xticklabels(highsympointlabels, fontsize = 15)
        ax[0].set_ylabel("Energy (eV)", fontsize = 20)
    
    
        ###DOS plot
        # Plot the wannier projected density of states
        xdos, TMdos, Liganddos = readwandos(dft_path, TM, TMorb, Ligand, Ligandorb)
    
        xmax = max(np.max(Liganddos), np.max(TMdos)) * 1.10
    
        ax[1].plot(TMdos,xdos,color=TMcolor)
        ax[1].plot(Liganddos,xdos,color=Ligandcolor)
    
        ax[1].set_ylim((emin,emax))
        ax[1].set_xlim((np.min(TMdos),xmax))
    
        # Set exactly 10 ticks, for aesthetics
        num_ticks = 10
        locator = ticker.MaxNLocator(num_ticks)
        ax[1].xaxis.set_major_locator(locator)
        ticks = ax[1].get_xticks()
        ax[1].xaxis.set_major_locator(ticker.FixedLocator(ticks))
    
        ax[1].set_xlabel("L DOS", labelpad = 4, fontsize = 15)
    
        ax[1].text(np.min(TMdos)+0.5*(xmax-np.min(TMdos)), emin+0.94*(emax-emin),
                   CF.Convert_CompName_ToLaTeX(compound_name), fontsize = 20, weight='bold',
                   horizontalalignment='center', verticalalignment='center')
    
        ax[0].tick_params(axis='both', which='both', direction='in', bottom=False, top=False, left=True, right=True,
               labelbottom=True, labeltop=False, labelleft=True, labelright=False)
        ax[1].tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True,
               labelbottom=False, labeltop=False, labelleft=False, labelright=False)
    
    
        TMdostemp = [x for x in TMdos if x>0.01]
        Liganddostemp = [x for x in TMdos if x>0.01]
    
        temp = TMtoplot
        TMavg = np.average(ywan[temp[0]:temp[-1]][0])
        temp = Ligandtoplot
        Ligandavg = np.average(ywan[temp[0]:temp[-1]][0])
    
        ax[1].text(np.average(TMdostemp)*4.5,TMavg+2,TM+"-"+TMorb,fontsize=12,color=TMcolor,style='italic')
        ax[1].text(np.average(Liganddostemp)*4.5,Ligandavg+2,Ligand+"-"+Ligandorb,fontsize=12,color=Ligandcolor,style='italic')
    
        plt.savefig(output_filename, dpi = 200, bbox_inches = 'tight')
        plt.close()

    makeBandStructDOSplot(dft_path = Path(dft_path), 
                          output_filename = plot_params_dict["output_filename"], 
                          emin = plot_params_dict["emin"],
                          emax = plot_params_dict["emax"], 
                          compound_name = material_input["NAME"],
                          TM = material_input["TM"], 
                          Ligand = material_input["LIGAND"], 
                          TMorb = material_input["TM_VALENCE"],
                          Ligandorb = material_input["LIGAND_VALENCE"])



def make_QuantyPlots_XES(data_path, plot_params_dict):
    """
    Automatically create summary figure of CTC-XES results,
    assumes the all linear and circular polarizations have been calculated
    """

    fig, ax = plt.subplots(4, 3, sharex=True, figsize=(12,7), gridspec_kw={'height_ratios':[3,1,3,1]})
    fig.subplots_adjust(hspace=0, wspace=0)
    data_path = Path(data_path)
    
    LEFT = plot_params_dict["left"]
    RIGHT = plot_params_dict["right"]

    # Calculate x_shift
    data_x = extract_quanty_spectra(data_path / "XES_xpol.dat")
    data_y = extract_quanty_spectra(data_path / "XES_ypol.dat")
    data_z = extract_quanty_spectra(data_path / "XES_zpol.dat")
    data_iso = np.asarray([data_x[0], (data_x[1] + data_y[1] + data_z[1]) / 3])
    X_SHIFT = -1 * data_iso[0][np.argmax(data_iso[1])]
    
    # Read in and process quanty spectra
    data_z_lin = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_zpol.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_y_lin = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_ypol.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_x_lin = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_xpol.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    
    data_z_rig = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_rpol_zdir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_z_lef = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_lpol_zdir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_y_rig = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_rpol_ydir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_y_lef = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_lpol_ydir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_x_rig = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_rpol_xdir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_x_lef = process_quanty_spectra_xes(extract_quanty_spectra(data_path / "XES_lpol_xdir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    
    data_iso = process_quanty_spectra_xes(data_iso, x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    x_default = data_z_lin[0]
    
    ### Z Propegation dir
    # linear
    ax[0][0].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][0].plot(data_x_lin[0], data_x_lin[1], label = "(z,x,0)", **plot_params_dict["x_lin"])
    ax[0][0].plot(data_y_lin[0], data_y_lin[1], label = "(z,y,0)", **plot_params_dict["y_lin"])
    ax[0][0].legend(loc=2)
    
    linear_diff = data_x_lin[1] - data_y_lin[1] # Linear Dichroism
    ax[1][0].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "X - Y")
    ax[1][0].legend(loc=2, handlelength = 1.0)
    
    # circular
    ax[2][0].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[2][0].plot(data_z_rig[0], data_z_rig[1], label = "(z,0,+1)", **plot_params_dict["z_rig"])
    ax[2][0].plot(data_z_lef[0], data_z_lef[1], label = "(z,0,-1)", **plot_params_dict["z_lef"])
    ax[2][0].legend(loc=2)
    
    circular_diff = data_z_rig[1] - data_z_lef[1] # Circular Dichroism
    ax[3][0].fill_between(x_default, np.zeros(len(circular_diff)), circular_diff, label = "R - L")
    ax[3][0].legend(loc=2, handlelength = 1.0)
    
    ### Y Propegation dir
    # linear
    ax[0][1].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][1].plot(data_z_lin[0], data_z_lin[1], label = "(y,z,0)", **plot_params_dict["z_lin"])
    ax[0][1].plot(data_x_lin[0], data_x_lin[1], label = "(y,x,0)", **plot_params_dict["x_lin"])
    ax[0][1].legend(loc=2)
    
    linear_diff = data_z_lin[1] - data_x_lin[1] # Linear Dichroism
    ax[1][1].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "Z - X")
    ax[1][1].legend(loc=2, handlelength = 1.0)
    
    # circular
    ax[2][1].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[2][1].plot(data_y_rig[0], data_y_rig[1], label = "(y,0,+1)", **plot_params_dict["y_rig"])
    ax[2][1].plot(data_y_lef[0], data_y_lef[1], label = "(y,0,-1)", **plot_params_dict["y_lef"])
    ax[2][1].legend(loc=2)
    
    circular_diff = data_y_rig[1] - data_y_lef[1] # Circular Dichroism
    ax[3][1].fill_between(x_default, np.zeros(len(circular_diff)), circular_diff, label = "R - L")
    ax[3][1].legend(loc=2, handlelength = 1.0)
    
    ### X Propegation dir
    # linear
    ax[0][2].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][2].plot(data_y_lin[0], data_y_lin[1], label = "(x,y,0)", **plot_params_dict["y_lin"])
    ax[0][2].plot(data_z_lin[0], data_z_lin[1], label = "(x,z,0)", **plot_params_dict["z_lin"])
    ax[0][2].legend(loc=2)
    
    linear_diff = data_y_lin[1] - data_z_lin[1] # Linear Dichroism
    ax[1][2].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "Y - Z")
    ax[1][2].legend(loc=2, handlelength = 1.0)
    
    # circular
    ax[2][2].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[2][2].plot(data_x_rig[0], data_x_rig[1], label = "(x,0,+1)", **plot_params_dict["x_rig"])
    ax[2][2].plot(data_x_lef[0], data_x_lef[1], label = "(x,0,-1)", **plot_params_dict["x_lef"])
    ax[2][2].legend(loc=2)
    
    circular_diff = data_x_rig[1] - data_x_lef[1] # Circular Dichroism
    ax[3][2].fill_between(x_default, np.zeros(len(circular_diff)), circular_diff, label = "R - L")
    ax[3][2].legend(loc=2, handlelength = 1.0)
    
    ### Aesthetics
    max_y = -float('inf')
    for subplot in ax[0]:
        ylim = subplot.get_ylim()
        max_y = max(max_y, ylim[1])
    for subplot in ax[2]:
        ylim = subplot.get_ylim()
        max_y = max(max_y, ylim[1])

    ax[0][1].set_yticklabels([])
    ax[0][2].set_yticklabels([])
    ax[0][1].sharey(ax[0][0])
    ax[0][2].sharey(ax[0][0])
    ax[0][0].set_ylim(None, max_y * 1.03)
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][0].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][0].get_yticks()]
    ax[1][0].set_yticks(ax[1][0].get_yticks()[1:-1])
    ax[1][0].set_yticklabels(tick_labels[1:-1])
    ax[1][0].text(0.05, 0.10, sn_label, transform=ax[1][0].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][1].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][1].get_yticks()]
    ax[1][1].set_yticks(ax[1][1].get_yticks()[1:-1])
    ax[1][1].set_yticklabels(tick_labels[1:-1])
    ax[1][1].text(0.05, 0.10, sn_label, transform=ax[1][1].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][2].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][2].get_yticks()]
    ax[1][2].set_yticks(ax[1][2].get_yticks()[1:-1])
    ax[1][2].set_yticklabels(tick_labels[1:-1])
    ax[1][2].text(0.05, 0.10, sn_label, transform=ax[1][2].transAxes, ha='left', va='bottom')
    
    ax[2][1].set_yticklabels([])
    ax[2][2].set_yticklabels([])
    ax[2][1].sharey(ax[2][0])
    ax[2][2].sharey(ax[2][0])
    ax[2][0].set_ylim(None, max_y * 1.03)
    
    exponent = np.abs(np.floor(np.log10(max(ax[3][0].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[3][0].get_yticks()]
    ax[3][0].set_yticks(ax[3][0].get_yticks()[1:-1])
    ax[3][0].set_yticklabels(tick_labels[1:-1])
    ax[3][0].text(0.05, 0.05, sn_label, transform=ax[3][0].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[3][1].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[3][1].get_yticks()]
    ax[3][1].set_yticks(ax[3][1].get_yticks()[1:-1])
    ax[3][1].set_yticklabels(tick_labels[1:-1])
    ax[3][1].text(0.05, 0.05, sn_label, transform=ax[3][1].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[3][2].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[3][2].get_yticks()]
    ax[3][2].set_yticks(ax[3][2].get_yticks()[1:-1])
    ax[3][2].set_yticklabels(tick_labels[1:-1])
    ax[3][2].text(0.05, 0.05, sn_label, transform=ax[3][2].transAxes, ha='left', va='bottom')

    # Axes Labels, Font Size, and Title
    title = CF.Convert_CompName_ToLaTeX(plot_params_dict["NAME"])
    if "Ka" in plot_params_dict["output_filename"]:
        fig.supxlabel(r'E - E$_{\mathdefault{K \alpha 1}}$ (eV)', fontsize = 20, y = 0.02)
        title += r" K$\alpha$"
    if "Kb" in plot_params_dict["output_filename"]:
        fig.supxlabel(r'E - E$_{\mathdefault{K \beta 1}}$ (eV)', fontsize = 20, y = 0.02)
        title += r" K$\beta$"

    title += " Polarized CTC-XES"
    fig.suptitle(title, fontsize = 30, y = 0.95)
    fig.supylabel('Intensity (arb. units)', fontsize = 20, x = 0.05)
    
    plt.savefig(plot_params_dict["output_filename"], dpi=200, bbox_inches='tight')
    plt.close()

def make_QuantyPlots_XAS(data_path, plot_params_dict):
    """
    Automatically create summary figure of XAS results,
    assumes the all linear and circular polarizations have been calculated
    """

    fig, ax = plt.subplots(4, 3, sharex=True, figsize=(12,7), gridspec_kw={'height_ratios':[3,1,3,1]})
    fig.subplots_adjust(hspace=0, wspace=0)
    data_path = Path(data_path)
    
    LEFT = plot_params_dict["left"]
    RIGHT = plot_params_dict["right"]

    # Calculate x_shift
    data_x = extract_quanty_spectra(data_path / "XAS_xpol.dat")
    data_y = extract_quanty_spectra(data_path / "XAS_ypol.dat")
    data_z = extract_quanty_spectra(data_path / "XAS_zpol.dat")
    data_iso = np.asarray([data_x[0], (data_x[1] + data_y[1] + data_z[1]) / 3])
    X_SHIFT = -1 * data_iso[0][np.argmax(data_iso[1])]
    
    # Read in and process quanty spectra
    data_z_lin = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_zpol.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_y_lin = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_ypol.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_x_lin = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_xpol.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    
    data_z_rig = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_rpol_zdir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_z_lef = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_lpol_zdir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_y_rig = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_rpol_ydir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_y_lef = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_lpol_ydir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_x_rig = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_rpol_xdir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_x_lef = process_quanty_spectra_xas(extract_quanty_spectra(data_path / "XAS_lpol_xdir.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    
    data_iso = process_quanty_spectra_xas(data_iso, x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    x_default = data_z_lin[0]
    
    ### Z Propegation dir
    # linear
    ax[0][0].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][0].plot(data_x_lin[0], data_x_lin[1], label = "(z,x,0)", **plot_params_dict["x_lin"])
    ax[0][0].plot(data_y_lin[0], data_y_lin[1], label = "(z,y,0)", **plot_params_dict["y_lin"])
    ax[0][0].legend(loc=2)
    
    linear_diff = data_x_lin[1] - data_y_lin[1] # Linear Dichroism
    ax[1][0].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "X - Y")
    ax[1][0].legend(loc=2, handlelength = 1.0)
    
    # circular
    ax[2][0].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[2][0].plot(data_z_rig[0], data_z_rig[1], label = "(z,0,+1)", **plot_params_dict["z_rig"])
    ax[2][0].plot(data_z_lef[0], data_z_lef[1], label = "(z,0,-1)", **plot_params_dict["z_lef"])
    ax[2][0].legend(loc=2)
    
    circular_diff = data_z_rig[1] - data_z_lef[1] # Circular Dichroism
    ax[3][0].fill_between(x_default, np.zeros(len(circular_diff)), circular_diff, label = "R - L")
    ax[3][0].legend(loc=2, handlelength = 1.0)
    
    ### Y Propegation dir
    # linear
    ax[0][1].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][1].plot(data_z_lin[0], data_z_lin[1], label = "(y,z,0)", **plot_params_dict["z_lin"])
    ax[0][1].plot(data_x_lin[0], data_x_lin[1], label = "(y,x,0)", **plot_params_dict["x_lin"])
    ax[0][1].legend(loc=2)
    
    linear_diff = data_z_lin[1] - data_x_lin[1] # Linear Dichroism
    ax[1][1].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "Z - X")
    ax[1][1].legend(loc=2, handlelength = 1.0)
    
    # circular
    ax[2][1].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[2][1].plot(data_y_rig[0], data_y_rig[1], label = "(y,0,+1)", **plot_params_dict["y_rig"])
    ax[2][1].plot(data_y_lef[0], data_y_lef[1], label = "(y,0,-1)", **plot_params_dict["y_lef"])
    ax[2][1].legend(loc=2)
    
    circular_diff = data_y_rig[1] - data_y_lef[1] # Circular Dichroism
    ax[3][1].fill_between(x_default, np.zeros(len(circular_diff)), circular_diff, label = "R - L")
    ax[3][1].legend(loc=2, handlelength = 1.0)
    
    ### X Propegation dir
    # linear
    ax[0][2].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][2].plot(data_y_lin[0], data_y_lin[1], label = "(x,y,0)", **plot_params_dict["y_lin"])
    ax[0][2].plot(data_z_lin[0], data_z_lin[1], label = "(x,z,0)", **plot_params_dict["z_lin"])
    ax[0][2].legend(loc=2)
    
    linear_diff = data_y_lin[1] - data_z_lin[1] # Linear Dichroism
    ax[1][2].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "Y - Z")
    ax[1][2].legend(loc=2, handlelength = 1.0)
    
    # circular
    ax[2][2].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[2][2].plot(data_x_rig[0], data_x_rig[1], label = "(x,0,+1)", **plot_params_dict["x_rig"])
    ax[2][2].plot(data_x_lef[0], data_x_lef[1], label = "(x,0,-1)", **plot_params_dict["x_lef"])
    ax[2][2].legend(loc=2)
    
    circular_diff = data_x_rig[1] - data_x_lef[1] # Circular Dichroism
    ax[3][2].fill_between(x_default, np.zeros(len(circular_diff)), circular_diff, label = "R - L")
    ax[3][2].legend(loc=2, handlelength = 1.0)
    
    ### Aesthetics
    max_y = -float('inf')
    for subplot in ax[0]:
        ylim = subplot.get_ylim()
        max_y = max(max_y, ylim[1])
    for subplot in ax[2]:
        ylim = subplot.get_ylim()
        max_y = max(max_y, ylim[1])

    ax[0][1].set_yticklabels([])
    ax[0][2].set_yticklabels([])
    ax[0][1].sharey(ax[0][0])
    ax[0][2].sharey(ax[0][0])
    ax[0][0].set_ylim(None, max_y * 1.03)
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][0].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][0].get_yticks()]
    ax[1][0].set_yticks(ax[1][0].get_yticks()[1:-1])
    ax[1][0].set_yticklabels(tick_labels[1:-1])
    ax[1][0].text(0.05, 0.10, sn_label, transform=ax[1][0].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][1].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][1].get_yticks()]
    ax[1][1].set_yticks(ax[1][1].get_yticks()[1:-1])
    ax[1][1].set_yticklabels(tick_labels[1:-1])
    ax[1][1].text(0.05, 0.10, sn_label, transform=ax[1][1].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][2].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][2].get_yticks()]
    ax[1][2].set_yticks(ax[1][2].get_yticks()[1:-1])
    ax[1][2].set_yticklabels(tick_labels[1:-1])
    ax[1][2].text(0.05, 0.10, sn_label, transform=ax[1][2].transAxes, ha='left', va='bottom')
    
    ax[2][1].set_yticklabels([])
    ax[2][2].set_yticklabels([])
    ax[2][1].sharey(ax[2][0])
    ax[2][2].sharey(ax[2][0])
    ax[2][0].set_ylim(None, max_y * 1.03)
    
    exponent = np.abs(np.floor(np.log10(max(ax[3][0].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[3][0].get_yticks()]
    ax[3][0].set_yticks(ax[3][0].get_yticks()[1:-1])
    ax[3][0].set_yticklabels(tick_labels[1:-1])
    ax[3][0].text(0.05, 0.05, sn_label, transform=ax[3][0].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[3][1].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[3][1].get_yticks()]
    ax[3][1].set_yticks(ax[3][1].get_yticks()[1:-1])
    ax[3][1].set_yticklabels(tick_labels[1:-1])
    ax[3][1].text(0.05, 0.05, sn_label, transform=ax[3][1].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[3][2].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[3][2].get_yticks()]
    ax[3][2].set_yticks(ax[3][2].get_yticks()[1:-1])
    ax[3][2].set_yticklabels(tick_labels[1:-1])
    ax[3][2].text(0.05, 0.05, sn_label, transform=ax[3][2].transAxes, ha='left', va='bottom')

    # Axes Labels, Font Size, and Title
    title = CF.Convert_CompName_ToLaTeX(plot_params_dict["NAME"])
    if "L23" in plot_params_dict["output_filename"]:
        title += r" L23 XAS"

    fig.supxlabel('Energy (eV)', fontsize = 20, y = 0.02)
    title += " Polarized XAS"
    fig.suptitle(title, fontsize = 30, y = 0.95)
    fig.supylabel('Intensity (arb. units)', fontsize = 20, x = 0.05)
    
    plt.savefig(plot_params_dict["output_filename"], dpi=200, bbox_inches='tight')
    plt.close()




def make_FEFFPlots(data_path, plot_params_dict):
    """
    Automatically create summary figure of VTC-XES results,
    assumes that isotropic and all linear polarizations have
    been calculated
    """

    fig, ax = plt.subplots(3, 3, sharex=True, figsize=(12,7), gridspec_kw={'height_ratios':[3,1,3]})
    fig.subplots_adjust(hspace=0, wspace=0)
    data_path = Path(data_path)
    
    LEFT = plot_params_dict["left"]
    RIGHT = plot_params_dict["right"]

    # Calculate x_shift
    data_iso = extract_feff_spectra(data_path / "isotropic" / "xmu.dat")
    X_SHIFT = -1 * data_iso[0][np.argmax(data_iso[1])]
    
    ### Read in and process feff spectra
    # Just Dipole
    data_z_lin_jd = process_feff_spectra(extract_feff_spectra(data_path / "z_polarization" / "xmu.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_y_lin_jd = process_feff_spectra(extract_feff_spectra(data_path / "y_polarization" / "xmu.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_x_lin_jd = process_feff_spectra(extract_feff_spectra(data_path / "x_polarization" / "xmu.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)

    data_x_dir_jd = avg_two_spectra(data_y_lin_jd, data_z_lin_jd)
    data_y_dir_jd = avg_two_spectra(data_z_lin_jd, data_x_lin_jd)
    data_z_dir_jd = avg_two_spectra(data_x_lin_jd, data_y_lin_jd)

    data_iso_jd = process_feff_spectra(data_iso, x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    x_default = data_z_lin_jd[0] # Energy values (x-axis), NOT x spatially (direction or polarization)

    # Dipole + Quadrupole
    data_iso_dq = process_feff_spectra(extract_feff_spectra(data_path / "isotropic_withquad" / "xmu.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_z_dir_dq = process_feff_spectra(extract_feff_spectra(data_path / "z_direction_withquad" / "xmu.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_y_dir_dq = process_feff_spectra(extract_feff_spectra(data_path / "y_direction_withquad" / "xmu.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)
    data_x_dir_dq = process_feff_spectra(extract_feff_spectra(data_path / "x_direction_withquad" / "xmu.dat"), x_shift = X_SHIFT, left = LEFT, right = RIGHT)

    # Just Quadrupole
    data_iso_jq = np.asarray([x_default, data_iso_dq[1] - data_iso_jd[1]])
    data_x_dir_jq = np.asarray([x_default, data_x_dir_dq[1] - data_x_dir_jd[1]])
    data_y_dir_jq = np.asarray([x_default, data_y_dir_dq[1] - data_y_dir_jd[1]])
    data_z_dir_jq = np.asarray([x_default, data_z_dir_dq[1] - data_z_dir_jd[1]])

   
    ### Z Propegation dir
    # linear
    ax[0][0].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][0].plot(data_x_lin_jd[0], data_x_lin_jd[1], label = "(z,x,0)", **plot_params_dict["x_lin"])
    ax[0][0].plot(data_y_lin_jd[0], data_y_lin_jd[1], label = "(z,y,0)", **plot_params_dict["y_lin"])
    ax[0][0].legend(loc=2)
    
    linear_diff = data_x_lin_jd[1] - data_y_lin_jd[1] # Linear Dichroism
    ax[1][0].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "X - Y")
    ax[1][0].legend(loc=2, handlelength = 1.0)

    ax[2][0].fill_between(data_z_dir_jd[0], np.zeros(len(data_z_dir_jd[1])), data_z_dir_jd[1], label = "[z,d,0]", **plot_params_dict["z_dir_jd"])
    ax[2][0].fill_between(data_z_dir_jq[0], np.zeros(len(data_z_dir_jq[1])), data_z_dir_jq[1], label = "[z,0,q]", **plot_params_dict["z_dir_jq"])
    ax[2][0].plot(data_z_dir_dq[0], data_z_dir_dq[1], label = "[z,d,q]", **plot_params_dict["z_dir_dq"])
    ax[2][0].legend(loc=2)
    
    ### Y Propegation dir
    # linear
    ax[0][1].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][1].plot(data_z_lin_jd[0], data_z_lin_jd[1], label = "(y,z,0)", **plot_params_dict["z_lin"])
    ax[0][1].plot(data_x_lin_jd[0], data_x_lin_jd[1], label = "(y,x,0)", **plot_params_dict["x_lin"])
    ax[0][1].legend(loc=2)
    
    linear_diff = data_z_lin_jd[1] - data_x_lin_jd[1] # Linear Dichroism
    ax[1][1].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "Z - X")
    ax[1][1].legend(loc=2, handlelength = 1.0)

    ax[2][1].fill_between(data_y_dir_jd[0], np.zeros(len(data_y_dir_jd[1])), data_y_dir_jd[1], label = "[y,d,0]", **plot_params_dict["y_dir_jd"])
    ax[2][1].fill_between(data_y_dir_jq[0], np.zeros(len(data_y_dir_jq[1])), data_y_dir_jq[1], label = "[y,0,q]", **plot_params_dict["y_dir_jq"])
    ax[2][1].plot(data_y_dir_dq[0], data_y_dir_dq[1], label = "[y,d,q]", **plot_params_dict["y_dir_dq"])
    ax[2][1].legend(loc=2)
    
    ### X Propegation dir
    # linear
    ax[0][2].plot(data_iso[0], data_iso[1], label = "Iso", **plot_params_dict["iso"])
    ax[0][2].plot(data_y_lin_jd[0], data_y_lin_jd[1], label = "(x,y,0)", **plot_params_dict["y_lin"])
    ax[0][2].plot(data_z_lin_jd[0], data_z_lin_jd[1], label = "(x,z,0)", **plot_params_dict["z_lin"])
    ax[0][2].legend(loc=2)
    
    linear_diff = data_y_lin_jd[1] - data_z_lin_jd[1] # Linear Dichroism
    ax[1][2].fill_between(x_default, np.zeros(len(linear_diff)), linear_diff, label = "Y - Z")
    ax[1][2].legend(loc=2, handlelength = 1.0)
    
    ax[2][2].fill_between(data_x_dir_jd[0], np.zeros(len(data_x_dir_jd[1])), data_x_dir_jd[1], label = "[x,d,0]", **plot_params_dict["x_dir_jd"])
    ax[2][2].fill_between(data_x_dir_jq[0], np.zeros(len(data_x_dir_jq[1])), data_x_dir_jq[1], label = "[x,0,q]", **plot_params_dict["x_dir_jq"])
    ax[2][2].plot(data_x_dir_dq[0], data_x_dir_dq[1], label = "[x,d,q]", **plot_params_dict["x_dir_dq"])
    ax[2][2].legend(loc=2)

    ### Aesthetics
    max_y = -float('inf')
    for subplot in ax[0]:
        ylim = subplot.get_ylim()
        max_y = max(max_y, ylim[1])

    ax[0][1].set_yticklabels([])
    ax[0][2].set_yticklabels([])
    ax[0][1].sharey(ax[0][0])
    ax[0][2].sharey(ax[0][0])
    ax[0][0].set_ylim(None, max_y * 1.03)
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][0].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][0].get_yticks()]
    ax[1][0].set_yticks(ax[1][0].get_yticks()[1:-1])
    ax[1][0].set_yticklabels(tick_labels[1:-1])
    ax[1][0].text(0.05, 0.05, sn_label, transform=ax[1][0].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][1].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][1].get_yticks()]
    ax[1][1].set_yticks(ax[1][1].get_yticks()[1:-1])
    ax[1][1].set_yticklabels(tick_labels[1:-1])
    ax[1][1].text(0.05, 0.05, sn_label, transform=ax[1][1].transAxes, ha='left', va='bottom')
    
    exponent = np.abs(np.floor(np.log10(max(ax[1][2].get_yticks()))))
    sn_label = "x$10^{" + str(-1 * int(exponent)) + "}$"
    sn_label = r"{}".format(sn_label)
    tick_labels = ['{:.2f}'.format(val * 10**(exponent)) for val in ax[1][2].get_yticks()]
    ax[1][2].set_yticks(ax[1][2].get_yticks()[1:-1])
    ax[1][2].set_yticklabels(tick_labels[1:-1])
    ax[1][2].text(0.05, 0.05, sn_label, transform=ax[1][2].transAxes, ha='left', va='bottom')
   
    # Axes Labels, Font Size, and Title
    title = CF.Convert_CompName_ToLaTeX(plot_params_dict["NAME"])
    fig.supxlabel(r'Relative Energy (eV)', fontsize = 20, y = 0.02)
    title += " Polarized VTC-XES"
    fig.suptitle(title, fontsize = 30, y = 0.95)
    fig.supylabel('Intensity (arb. units)', fontsize = 20, x = 0.05)
 
    plt.savefig(plot_params_dict["output_filename"], dpi=200, bbox_inches='tight')
    plt.close()


def make_crystalvisplots(structure, plot_params_dict):
    """
    Generate 3 snapshots of the crystal structure using PyMatGen.
    Will create 3 images, along the x, y, and z axes.
    """

    # Tolerance to use when determining if atoms are at the same position
    TOL = 1e-5 
        
    frac_coords_list = [tuple(site.frac_coords) for site in structure.sites]
    sg_analyzer = SpacegroupAnalyzer(structure)
    symmetry_operations = sg_analyzer.get_symmetry_operations()
    
    new_structure = Structure.from_spacegroup(structure.get_space_group_info()[1], 
                                                  lattice=structure.lattice.matrix, 
                                                  species=[], coords=[])

    direction = plot_params_dict["direction"]
    
    def angle_between_vectors(v, w):
        dot_product = np.dot(v, w)
        norm_v = np.linalg.norm(v)
        norm_w = np.linalg.norm(w)
    
        # Ensure the denominators are not zero
        if norm_v == 0 or norm_w == 0:
            raise ValueError("Input vectors must have non-zero norms")
    
        cos_theta = dot_product / (norm_v * norm_w)
        angle_rad = np.arccos(np.clip(cos_theta, -1.0, 1.0))
        angle_deg = np.degrees(angle_rad)
    
        return angle_deg
    
    
    def are_vectors_close(vec1, vec2):
        if len(vec1) != len(vec2):
            raise ValueError("Vectors must have the same length")
    
        for elem1, elem2 in zip(vec1, vec2):
            if abs(elem1 - elem2) >= TOL:
                return False
        
        return True
    
    def site_exists(structure, species, position):
        for site in structure:
            coords = site.frac_coords
            species_string = site.specie.symbol
            if species == species_string and are_vectors_close(list(coords), position):
                return True
        return False
    
    
    # Apply a, b, c shifts to every atom in the structure
    for vec in [[0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,1,0], [0,1,1], [1,0,1], [1,1,1]]:
        for site in structure:
            new_position = np.asarray(vec) + site.frac_coords
            for key, val in enumerate(new_position):
                if np.abs(val - 1) < TOL:
                    new_position[key] = 1
                elif np.abs(val) < TOL:
                    new_position[key] = 0
            new_position = tuple(new_position)
            outside_uc = False
            for key, val in enumerate(new_position):
                if val < 0 or val > 1:
                    outside_uc = True
    
            if site_exists(new_structure, site.specie.symbol, list(new_position)) or outside_uc:
                continue
            else:
                new_structure.append(site.species_string, list(new_position))
    

    ### Make screenshot of structure when facing down a particular axes

    if direction not in ["x", "y", "z"]:
        raise ValueError(f"Invalid direction: {direction}")

    vis = StructureVis(show_bonds=True, show_polyhedron=False, show_unit_cell = True)
    vis.set_structure(new_structure, to_unit_cell=False, reset_camera=True)
    
    # Get a list of every atom coordinate that will be visualized
    all_atom_coords = []
    for key, site in vis.mapper_map.items():
        all_atom_coords.append(site[0].coords)
    
    # Manage Camera Position
    matrix = new_structure.lattice.matrix
    theta = angle_between_vectors(np.asarray([1,0,0]), matrix[0])
    camera = vis.ren.GetActiveCamera()
    
    # Adjust the camera for best viewing
    lengths = new_structure.lattice.abc
    pos = (matrix[1] + matrix[2]) * 0.5 + matrix[0] * max(lengths) / lengths[0] * 5
    camera.SetPosition(pos)
    camera.SetViewUp(matrix[2])
    camera.SetFocalPoint((matrix[0] + matrix[1] + matrix[2]) * 0.5)
    
    a, b, c = new_structure.lattice.abc
    if c > a:
        camera.SetViewAngle(50 - c/a * 5)
    else:
        camera.SetViewAngle(50 + c/a * 5)
    
    
    if direction == "x":
        camera.Elevation(theta)
        vis.ren_win.Render()
        vis.ren_win.SetSize(600, 600)
        
        # Set compass
        min_y = float('inf')
        for coord in all_atom_coords:
            if coord[1] < min_y:
                min_y = coord[1]
        
        arrow_factor = 3 # Controls width and length
        matrix = new_structure.lattice.matrix
        origin = (np.linalg.norm(matrix[0])/2, min_y - arrow_factor - new_structure.lattice.abc[1] * 0.3, np.linalg.norm(matrix[2])/2)
        origin = np.asarray(origin)
        vis.add_compass(origin, arrow_factor, arrow_factor)
        
    elif direction == "y":
        camera.Elevation(theta)
        camera.Azimuth(90)
        camera.Roll(-90)
        vis.ren_win.Render()
        vis.ren_win.SetSize(600, 600)
        
        # Set compass
        min_z = float('inf')
        for coord in all_atom_coords:
            if coord[2] < min_z:
                min_z = coord[2]
        
        arrow_factor = 3 # Controls width and length
        matrix = new_structure.lattice.matrix
        origin = (np.linalg.norm(matrix[0])/2, np.linalg.norm(matrix[1])/2, min_z - arrow_factor - new_structure.lattice.abc[2] * 0.3)
        origin = np.asarray(origin)
        vis.add_compass(origin, arrow_factor, arrow_factor)
        
    else:
        camera.Elevation(theta)
        camera.Elevation(90)
        camera.Roll(90)
        vis.ren_win.Render()
        vis.ren_win.SetSize(600, 600)
        
        # Set compass
        min_x = float('inf')
        for coord in all_atom_coords:
            if coord[0] < min_x:
                min_x = coord[0]
        
        arrow_factor = 3 # Controls width and length
        matrix = new_structure.lattice.matrix
        origin = (min_x - arrow_factor - new_structure.lattice.abc[0] * 0.3, np.linalg.norm(matrix[1])/2, np.linalg.norm(matrix[2])/2)
        origin = np.asarray(origin)
        vis.add_compass(origin, arrow_factor, arrow_factor)
        
    vis.write_image(filename=plot_params_dict["output_filename"])

def stitch_images_horizontal(image_paths, output_path):
    """
    Combines images by appending them horizontally into a
    single image
    """
    
    images = [Image.open(path) for path in image_paths]
    
    # Get the width and height of the first image
    total_width = sum(img.width for img in images)
    max_height = max(img.height for img in images)

    # Create a new image with the computed dimensions
    result_image = Image.new('RGB', (total_width, max_height))

    # Paste each image into the result image horizontally
    current_width = 0
    for img in images:
        result_image.paste(img, (current_width, 0))
        current_width += img.width

    # Save the final image
    result_image.save(output_path)

def crop_images(filePaths):
    """
    Crops images to remove extra whitespace
    """
    padding = np.asarray([-5, -5, 5, 5])
    if not isinstance(filePaths, list):
        filePaths = [filePaths]
    for filePath in filePaths:
        image=Image.open(filePath)
        image.load()
        imageSize = image.size
    
        # remove alpha channel
        invert_im = image.convert("RGB")
    
        # invert image (so that white is 0)
        invert_im = ImageOps.invert(invert_im)
        imageBox = invert_im.getbbox()
        imageBox = tuple(np.asarray(imageBox)+padding)
    
        cropped=image.crop(imageBox)
        print(filePath, "Size:", imageSize, "New Size:", (imageBox[2] - imageBox[0], imageBox[3] - imageBox[1]))
        cropped.save(filePath)


def convert_images_to_pdf(image_paths, output_path):
    pdf_writer = PdfWriter()

    for image_path in image_paths:
        # Open each image and convert it to RGB
        img = Image.open(image_path).convert("RGB")

        # Create a temporary PDF file for each image
        temp_pdf_path = os.path.splitext(image_path)[0] + ".pdf"
        img.save(temp_pdf_path, "PDF")

        # Merge the temporary PDF into the final PDF
        with open(temp_pdf_path, "rb") as temp_pdf_file:
            pdf_writer.append(fileobj=temp_pdf_file)

        # Remove the temporary PDF file
        os.remove(temp_pdf_path)

    # Save the final PDF file
    with open(output_path, "wb") as output_file:
        pdf_writer.write(output_file)

    print(f"Images converted to PDF and saved to {output_path}")


def move_png_files(folder_name):
    """
    Find all .png files in current folder and move
    them into a new folder. This is for cleaning up the results
    of all the above image generation.
    """

    # Check if the folder already exists, create it if not
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    # Get a list of all files in the current directory
    files = os.listdir()

    # Filter only .png files
    png_files = [file for file in files if file.lower().endswith('.png')]

    # Move each .png file to the new folder
    for png_file in png_files:
        old_path = os.path.abspath(png_file)
        new_path = os.path.join(folder_name, png_file)

        # Use shutil.move to perform the file move
        shutil.move(old_path, new_path)

    print(f"PNG files moved to the {folder_name} folder.")


def text_to_image(text, output_path='header.png', font_size=20, image_size=(2080, 1342), text_color=(0, 0, 0), background_color=(255, 255, 255)):

    # Create a new image with a white background
    img = Image.new('RGB', image_size, background_color)
    draw = ImageDraw.Draw(img)

    # Use the default font
    font = ImageFont.load_default(size = font_size)

    # Get the bounding box of the text
    bbox = draw.textbbox((0, 0), text, font=font)

    # Calculate the text position to center it both horizontally and vertically
    x = (image_size[0] - (bbox[2] - bbox[0])) // 2
    y = (image_size[1] - (bbox[3] - bbox[1])) // 2


    # Write the text to the image
    draw.text((x, y), text, font=font, fill=text_color)
 
    # Save the image to a file
    img.save(output_path)

